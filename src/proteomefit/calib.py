"""예측 신뢰도 보정 - 교차검증 예측 -> 기울기 부스팅 확률표

신뢰도 (name)
  ess  : 필수 유전자일 확률               입력 = 분류 로짓            보정 전 = 분류 모델 확률
  sign : 실측이 예측과 같은 쪽에 있을 확률  입력 = (mu, log sd)         보정 전 = Phi(|mu| / sd)
  fit  : 실측이 mu +- sd 안에 있을 확률     입력 = (mu, log sd)         보정 전 = 0.683
분할 - 외부 폴드 = 시험 genus, 내부 폴드 j = 원 모델 내부 폴드 j 의 검증 genus (최적 반복 수 선택) · 시험 = 원 모델 내부 폴드 j 의 예측,
       최종 반복 수 = genus별 중앙값의 중앙값
"""
import glob, os
import joblib
import numpy as np, pandas as pd
from scipy.interpolate import RegularGridInterpolator
from scipy.special import logit
from scipy.stats import norm
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import log_loss
from .cv import with_stats
from .metrics import calib_metrics_by_organism, in_total, reliability_table
from .plots import ess_plot_title, prob_confusion_plot, reliability_plot


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


def fit_hgb(X, hit, n_iter, mono, params):
    """입력 [m, d] · 적중 [m] · 반복 수 · 단조 제약 (None | [1]) · 부스팅 설정 dict -> 학습 모델 (조기 종료 없음)"""
    return HistGradientBoostingClassifier(max_iter=n_iter, early_stopping=False, monotonic_cst=mono,
                                          random_state=0, **params).fit(X, hit)


def staged_tables(model, axes):
    """학습 모델 · 격자 축 목록 -> 반복 수별 확률표 [n_iter, *격자 모양]"""
    G = np.column_stack([a.ravel() for a in np.meshgrid(*axes, indexing='ij')])   # 격자 점 [n_grid, d]
    return np.stack([p[:, 1].reshape([len(a) for a in axes]) for p in model.staged_predict_proba(G)])


def best_iter(model, X, hit):
    """학습 모델 · 검증 입력 [m, d] · 적중 [m] -> 검증 로그 손실이 가장 작은 반복 수 (1부터)"""
    return int(np.argmin([log_loss(hit, p[:, 1], labels=[0, 1]) for p in model.staged_predict_proba(X)])) + 1


def lookup(table, axes, X):
    """확률표 · 격자 축 목록 · 입력 [n, d] -> 확률 [n] (선형 보간, 격자 밖은 가장자리 값)"""
    Xc = np.column_stack([np.clip(X[:, i], a[0], a[-1]) for i, a in enumerate(axes)])
    return RegularGridInterpolator(tuple(axes), table)(Xc)


def nested_cv(name, df, val, axes, mono, params, max_iter, bins, out_dir):
    """신뢰도 이름 · 예측 표 · {시험 genus: 검증 genus 3개} · 격자 축 · 단조 제약 · 부스팅 설정 · 반복 상한 · 확률 구간 · 결과 폴더
       -> (genus별 결과 [n_fold, *], 시험 예측의 확률 구간별 표). 내부 폴드 j 시험 = genus 의 원 모델 내부 폴드 j 예측,
       genus 성능 = 내부 폴드 지표의 평균 · 표준편차, 내부 폴드 지표 = 균주별 지표의 평균 (ess 는 양 · 음성이 모두 있는 균주만)
       산출물 (out_dir 아래)
         fold##_<genus>/inner##/ : model.pkl (best_iter 회 학습 모델) · pred.npz (best_iter · 확률표 table · 시험 예측 after) ·
                                   metrics_organ.csv · calibration.png · confusion.png (ess). model.pkl · pred.npz 가 있으면 학습 없이 읽는다
         fold##_<genus>/         : metrics_folds.csv - 내부 폴드별 1줄 + mean · median · std 줄
         ./                      : metrics_total.csv - genus별 1줄 + genus 간 mean · median · std 줄"""
    X, hit, before = calib_target(name, df)
    gen, org, src = df.genus.values, df.org.values, df.inner.values
    rows, pool = [], []
    for k, (G, vs) in enumerate(val.items()):
        fdir, inn = f'{out_dir}/fold{k + 1:02d}_{G}', []
        for j, v in enumerate(vs):
            tr, va, te = (gen != G) & (gen != v), gen == v, (gen == G) & (src == j + 1)
            jd, tag = f'{fdir}/inner{j + 1:02d}', f'conf_{name} fold{k + 1:02d} {G} inner{j + 1:02d} (val {v})'
            if not os.path.exists(f'{jd}/pred.npz'):
                model = fit_hgb(X[tr], hit[tr], max_iter, mono, params)
                b     = best_iter(model, X[va], hit[va])               # 검증 genus 로 고른 반복 수
                table = staged_tables(model, axes)[b - 1]              # [*격자 모양]
                os.makedirs(jd, exist_ok=True)
                np.savez(f'{jd}/pred.npz', best_iter=b, table=table, after=lookup(table, axes, X[te]))
            z = np.load(f'{jd}/pred.npz')
            if not os.path.exists(f'{jd}/model.pkl'):                  # best_iter 회 학습 모델 (300 회 모델의 best_iter 단계와 같다)
                joblib.dump(fit_hgb(X[tr], hit[tr], int(z['best_iter']), mono, params), f'{jd}/model.pkl')
            b, probs = int(z['best_iter']), {'before': before[te], 'after': z['after']}   # 시험 확률 [n_te]
            om = calib_metrics_by_organism(org[te], hit[te], probs, bins, call=name == 'ess')
            om.round(6).to_csv(f'{jd}/metrics_organ.csv', index=False)
            om = om[in_total(om)] if name == 'ess' else om
            reliability_plot({f'conf_{name}': reliability_table(hit[te], probs, bins)}, f'{jd}/calibration.png', tag)
            if name == 'ess':
                after_cols = om.rename(columns=lambda c: c[:-6] if c.endswith('_after') else c)   # 보정 후 지표 -> 원래 이름
                prob_confusion_plot(hit[te], probs['after'], f'{jd}/confusion.png', ess_plot_title(f'{tag} | calibrated', int(te.sum()), after_cols),
                                    ('non-essential', 'essential'), 0.5)
            inn.append(dict(fold=k + 1, genus=G, inner=j + 1, val_genus=v, n_test=int(te.sum()), best_iter=b,
                            **om.drop(columns=['org_id', 'n', 'n_pos']).mean()))
            pool.append(pd.DataFrame({'before': probs['before'], 'after': probs['after'], 'hit': hit[te]}))
        inn = pd.DataFrame(inn)
        with_stats(inn, {'fold': k + 1, 'genus': G}, 'inner')[0].round(6).to_csv(f'{fdir}/metrics_folds.csv', index=False)
        mets = inn.columns.drop(['fold', 'genus', 'inner', 'val_genus', 'n_test', 'best_iter'])
        rows.append(dict(fold=k + 1, genus=G, val_genera=';'.join(vs), n_test=int(inn.n_test.mean()),
                         best_iter_med=inn.best_iter.median(), best_iter_mean=inn.best_iter.mean(), best_iter_sd=inn.best_iter.std(),
                         **{kv: x for c in mets for kv, x in ((c, inn[c].mean()), (c + '_sd', inn[c].std()))}))
        r = rows[-1]
        print(f"conf_{name} fold {k + 1:2d} {G:<16} best@{r['best_iter_med']:5.1f} | "
              f"ece {r['ece_before']:.3f} -> {r['ece_after']:.3f} | brier {r['brier_before']:.3f} -> {r['brier_after']:.3f}")
    res = pd.DataFrame(rows)
    with_stats(res, {'genus': ''}, 'fold')[0].round(6).to_csv(f'{out_dir}/metrics_total.csv', index=False)
    P = pd.concat(pool, ignore_index=True)
    return res, reliability_table(P.hit.values, {'before': P.before.values, 'after': P.after.values}, bins)


def final_model(name, df, axes, mono, params, n_iter, out_dir):
    """신뢰도 이름 · 예측 표 (전체 genus) · 격자 축 · 단조 제약 · 부스팅 설정 · 반복 수 · 결과 폴더
       -> 최종 확률표 [*격자 모양] (out_dir/model_final.pkl 저장)"""
    X, hit, _ = calib_target(name, df)
    model = fit_hgb(X, hit, n_iter, mono, params)
    joblib.dump(model, f'{out_dir}/model_final.pkl')
    return staged_tables(model, axes)[-1]
