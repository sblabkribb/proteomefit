"""예측 신뢰도 보정 - 교차검증 예측 -> 기울기 부스팅 확률표

신뢰도 (name)
  ess  : 필수 유전자일 확률               입력 = 분류 로짓            보정 전 = 분류 모델 확률
  sign : 실측이 예측과 같은 쪽에 있을 확률  입력 = (mu, log sd)         보정 전 = Phi(|mu| / sd)
  fit  : 실측이 mu +- sd 안에 있을 확률     입력 = (mu, log sd)         보정 전 = 0.683
분할 = 원 모델과 같다 - 외부 폴드 = 시험 genus, 내부 폴드 j = 원 모델 내부 폴드 j 의 검증 genus (최적 반복 수 best_epoch 선택),
       시험 = 원 모델 내부 폴드 j 의 예측, 최종 반복 수 = genus별 중앙값의 중앙값
산출물 (out_dir 아래) - 표의 열 이름 · 순서는 원 모델 표(cv.run_cv)와 같고, 보정 지표(ece · brier · log_loss)가 뒤에 붙는다
  fold##_<genus>/inner##/ : model_calib.pkl · pred_calib.npz (best_epoch · 확률표 table · 시험 확률 prob, 있으면 학습 없이 읽는다) ·
                            metrics_organ_calib.csv · calibration.png · confusion.png (ess, 보정 전 · 후)
  fold##_<genus>/         : metrics_folds_calib.csv - 내부 폴드별 1줄 + mean · median · std 줄
  ./                      : metrics_total_calib.csv - genus별 1줄 + genus 간 mean · median · std 줄
                            compare_calib.csv (보정 전 · 후 genus 평균 · 표준편차) · check_calib.csv · calibration.png · model_final_calib.pkl
"""
import glob, os
import joblib
import numpy as np, pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.special import logit
from scipy.stats import norm
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import log_loss
from .cv import genus_row, with_stats
from .metrics import CAL_COLS, ESS_METRIC_COLS, calib_metrics_by_organism, in_total, org_summary, reliability_table
from .plots import confusion_plot, reliability_plot

GRID_MU, GRID_LSD = np.round(np.arange(-8, 2.01, 0.1), 2), np.round(np.arange(-1.6, 1.21, 0.1), 2)
AXES     = {'ess': [np.round(np.arange(-10, 10.01, 0.1), 2)],                # 로짓 격자 [201]
            'sign': [GRID_MU, GRID_LSD], 'fit': [GRID_MU, GRID_LSD]}         # mu [101] · log(sd) [29]
MONO     = {'ess': [1], 'sign': None, 'fit': None}                          # 단조 제약 (ess = 로짓에 단조 증가)
SCORE    = {'ess': 'mcc', 'sign': 'log_loss', 'fit': 'log_loss'}            # 검증 지표 열 val_{SCORE} (균주 평균)
HGB      = dict(learning_rate=0.1, max_leaf_nodes=15, min_samples_leaf=1000)   # 기울기 부스팅 설정
MAX_ITER = 300                                                              # 내부 폴드 반복 상한
BINS     = np.linspace(0, 1, 11)                                            # 확률 구간 10개 (ECE · 확인 표)
ROUND    = 6                                                                # csv 실수 반올림 자릿수
CLASSES  = ('non-essential', 'essential')


def load_cv_preds(src, cols, gen):
    """교차검증 결과 폴더 · 샘플 열 dict {이름: [n]} · 샘플 genus [n] -> 예측 표 [n_pred, *] (pred.npz 의 열 + cols + genus · inner)
       pred.npz = 시험 genus 샘플 순서, 내부 폴드마다 1벌"""
    rows = []
    for d in sorted(glob.glob(f'{src}/fold*')):
        genus = os.path.basename(d).split('_', 1)[1]
        te    = np.flatnonzero(gen == genus)                           # 시험 샘플 [n_te]
        for j, f in enumerate(sorted(glob.glob(f'{d}/inner*/pred.npz'))):
            rows.append(pd.DataFrame({**np.load(f), **{c: v[te] for c, v in cols.items()}, 'genus': genus, 'inner': j + 1}))
    return pd.concat(rows, ignore_index=True)


def val_genera(src):
    """교차검증 결과 폴더 -> {시험 genus: [내부 검증 genus 3개]} (외부 폴드 순서, metrics_total.csv 의 val_genera)"""
    t = pd.read_csv(f'{src}/metrics_total.csv').dropna(subset=['genus'])
    return {g: v.split(';') for g, v in zip(t.genus, t.val_genera)}


def calib_target(name, df):
    """신뢰도 이름 · 예측 표 -> (입력 X [m, 1 | 2], 적중 hit [m] 0/1, 보정 전 확률 before [m])"""
    if name == 'ess':
        p = df.prob.values.astype(float)
        return logit(np.clip(p, 1e-6, 1 - 1e-6))[:, None], df.ess.values.astype(int), p
    mu, sd = df.mu.values.astype(float), df.sd.values.astype(float)
    X = np.column_stack([mu, np.log(sd)])
    if name == 'sign':
        return X, (np.sign(df.y) == np.sign(mu)).values.astype(int), norm.cdf(np.abs(mu) / sd)
    return X, (np.abs(df.y - mu) <= sd).values.astype(int), np.full(len(df), 0.683)


def fit_hgb(X, hit, n_iter, mono):
    """입력 [m, d] · 적중 [m] · 반복 수 · 단조 제약 (None | [1]) -> 학습 모델 (조기 종료 없음)"""
    return HistGradientBoostingClassifier(max_iter=n_iter, early_stopping=False, monotonic_cst=mono,
                                          random_state=0, **HGB).fit(X, hit)


def staged_tables(model, axes):
    """학습 모델 · 격자 축 목록 -> 반복 수별 확률표 [n_iter, *격자 모양]"""
    G = np.column_stack([a.ravel() for a in np.meshgrid(*axes, indexing='ij')])   # 격자 점 [n_grid, d]
    return np.stack([p[:, 1].reshape([len(a) for a in axes]) for p in model.staged_predict_proba(G)])


def best_epoch(model, X, hit):
    """학습 모델 · 검증 입력 [m, d] · 적중 [m] -> 검증 로그 손실이 가장 작은 반복 수 (1부터)"""
    return int(np.argmin([log_loss(hit, p[:, 1], labels=[0, 1]) for p in model.staged_predict_proba(X)])) + 1


def lookup(table, axes, X):
    """확률표 · 격자 축 목록 · 입력 [n, d] -> 확률 [n] (선형 보간, 격자 밖은 가장자리 값)"""
    Xc = np.column_stack([np.clip(X[:, i], a[0], a[-1]) for i, a in enumerate(axes)])
    return RegularGridInterpolator(tuple(axes), table)(Xc)


def org_table(name, org, hit, prob):
    """신뢰도 이름 · 균주 · 적중 · 확률 [m] -> (균주별 지표 [n_org, *], 성능 집계에 쓰는 균주 (ess = 양 · 음성이 모두 있는 균주))"""
    om = calib_metrics_by_organism(org, hit, prob, BINS, CAL_COLS[name])
    return om, (om[in_total(om)] if name == 'ess' else om)


def fit_inner(name, X, hit, tr, va, te, jd):
    """내부 폴드 학습 (pred_calib.npz · model_calib.pkl 이 있으면 읽기만) -> (best_epoch, 확률표, 검증 확률 [n_va], 시험 확률 [n_te])"""
    axes, mono = AXES[name], MONO[name]
    if not os.path.exists(f'{jd}/pred_calib.npz'):
        model = fit_hgb(X[tr], hit[tr], MAX_ITER, mono)
        b     = best_epoch(model, X[va], hit[va])                       # 검증 genus 로 고른 반복 수
        table = staged_tables(model, axes)[b - 1]                      # [*격자 모양]
        np.savez(f'{jd}/pred_calib.npz', best_epoch=b, table=table, prob=lookup(table, axes, X[te]))
    z = np.load(f'{jd}/pred_calib.npz')
    b = int(z['best_epoch'])
    if not os.path.exists(f'{jd}/model_calib.pkl'):                    # best_epoch 회 학습 모델 (MAX_ITER 회 모델의 best_epoch 단계와 같다)
        joblib.dump(fit_hgb(X[tr], hit[tr], b, mono), f'{jd}/model_calib.pkl')
    return b, z['table'], lookup(z['table'], axes, X[va]), z['prob']


def note(t):
    """균주별 지표 표 (또는 genus 표) -> 보정 그림 패널 지표 줄 (평균)"""
    return f'ECE {t.ece.mean():.3f}  Brier {t.brier.mean():.3f}  log loss {t.log_loss.mean():.3f}'


def nested_cv(name, df, src, out_dir):
    """신뢰도 이름 · 예측 표 · 원 모델 결과 폴더 · 결과 폴더 -> (genus별 보정 후 결과 [n_fold, *], 보정 전 · 후 비교 표, 확인 표)
       내부 폴드 행 = 원 모델 metrics_folds.csv 의 분할 열 (val_genus · val_sim · val_nn) + 보정 결과"""
    X, hit, before = calib_target(name, df)
    gen, org, src_inner = df.genus.values, df.org.values, df.inner.values
    cols, sc = CAL_COLS[name], f'val_{SCORE[name]}'
    base = ESS_METRIC_COLS if name == 'ess' else []                      # 원 모델 표와 같은 지표 (나머지 = 보정 지표, 뒤에 붙인다)
    rows, rows_b, pool = [], [], []
    for k, (G, vs) in enumerate(val_genera(src).items()):
        fdir = f'fold{k + 1:02d}_{G}'
        split = pd.read_csv(f'{src}/{fdir}/metrics_folds.csv').dropna(subset=['val_genus'])   # 원 모델 내부 폴드 행 [3, *]
        inn, inn_b = [], []
        for j, v in enumerate(vs):
            tr, va, te = (gen != G) & (gen != v), gen == v, (gen == G) & (src_inner == j + 1)
            jd, tag = f'{out_dir}/{fdir}/inner{j + 1:02d}', f'conf_{name} fold{k + 1:02d} {G} inner{j + 1:02d} (val {v})'
            os.makedirs(jd, exist_ok=True)
            b, _, p_va, p_te = fit_inner(name, X, hit, tr, va, te, jd)
            om, tot     = org_table(name, org[te], hit[te], p_te)            # 보정 후 균주별 지표
            _, tot_b    = org_table(name, org[te], hit[te], before[te])      # 보정 전 (비교 표용)
            om.insert(0, 'fold', k + 1); om.insert(1, 'genus', G); om.insert(2, 'inner', j + 1)
            om.round(ROUND).to_csv(f'{jd}/metrics_organ_calib.csv', index=False)
            s = split.iloc[j]
            inn.append(dict(fold=k + 1, genus=G, inner=j + 1, val_genus=v, val_sim=s.val_sim, val_nn=s.val_nn, n_val=int(va.sum()),
                            n_org=len(tot), best_epoch=b, **{sc: org_table(name, org[va], hit[va], p_va)[1][SCORE[name]].mean()},
                            **org_summary(tot, base), **org_summary(tot, cols[len(base):])))
            inn_b.append({c: tot_b[c].mean() for c in cols})
            reliability_plot({f'conf_{name}': reliability_table(hit[te], {'before': before[te], 'after': p_te}, BINS)}, f'{jd}/calibration.png', tag,
                             {f'conf_{name}': {'before': note(tot_b), 'after': note(tot)}}, BINS)
            if name == 'ess':
                confusion_plot(hit[te], {f'{st} calibration\nMCC {t.mcc.mean():.3f}  recall {t.recall.mean():.3f}  precision {t.precision.mean():.3f}':
                                         (p >= 0.5).astype(int) for st, p, t in (('before', before[te], tot_b), ('after', p_te, tot))},
                               f'{jd}/confusion.png', tag, CLASSES)
            pool.append(pd.DataFrame({'before': before[te], 'after': p_te, 'hit': hit[te]}))
        inn = pd.DataFrame(inn)
        with_stats(inn, {'fold': k + 1, 'genus': G}, 'inner')[0].round(ROUND).to_csv(f'{out_dir}/{fdir}/metrics_folds_calib.csv', index=False)
        extra = {'n_pos': int(hit[(gen == G) & (src_inner == 1)].sum())} if name == 'ess' else {}
        rows.append(genus_row(k, G, vs, int(((gen == G) & (src_inner == 1)).sum()), extra, inn, [sc] + cols))
        rows_b.append(pd.DataFrame(inn_b).mean())
        r = rows[-1]
        print(f"conf_{name} fold {k + 1:2d} {G:<16} best@{r['best_epoch_med']:5.1f} | {sc} {r[sc]:.3f} | "
              f"ece {rows_b[-1]['ece']:.3f} -> {r['ece']:.3f} | brier {rows_b[-1]['brier']:.3f} -> {r['brier']:.3f}")
    res, res_b = pd.DataFrame(rows), pd.DataFrame(rows_b)
    with_stats(res, {'genus': ''}, 'fold')[0].round(ROUND).to_csv(f'{out_dir}/metrics_total_calib.csv', index=False)
    cmp = pd.DataFrame([dict(metric=c, before=res_b[c].mean(), before_sd=res_b[c].std(), after=res[c].mean(), after_sd=res[c].std())
                        for c in cols])                                  # genus 평균 · 표준편차 (NaN genus 제외)
    cmp.round(4).to_csv(f'{out_dir}/compare_calib.csv', index=False)
    P = pd.concat(pool, ignore_index=True)
    check = reliability_table(P.hit.values, {'before': P.before.values, 'after': P.after.values}, BINS)
    check.round(4).to_csv(f'{out_dir}/check_calib.csv', index=False)
    reliability_plot({f'conf_{name}': check}, f'{out_dir}/calibration.png', f'conf_{name} - {os.path.basename(os.path.normpath(out_dir))} (genus mean)',
                     {f'conf_{name}': {'before': note(res_b), 'after': note(res)}}, BINS)
    return res, cmp, check


def final_model(name, df, res, out_dir):
    """신뢰도 이름 · 예측 표 (전체 genus) · nested_cv 결과 · 결과 폴더 -> 최종 확률표 dict {axes, p, n_iter} (model_final_calib.pkl 저장)
       반복 수 = genus별 best_epoch 중앙값의 중앙값"""
    X, hit, _ = calib_target(name, df)
    n_iter = int(res.best_epoch_med.median() + 0.5)
    model  = fit_hgb(X, hit, n_iter, MONO[name])
    joblib.dump(model, f'{out_dir}/model_final_calib.pkl')
    return dict(axes=[a.tolist() for a in AXES[name]], p=staged_tables(model, AXES[name])[-1].round(4).tolist(), n_iter=n_iter)
