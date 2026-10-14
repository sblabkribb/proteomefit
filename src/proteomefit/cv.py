"""genus 단위 중첩 교차검증 실행 · 최종 모델 · 결과 요약

산출물 (out_dir 아래)
  fold##_<genus>/inner##/ : task.artifacts (model.pt · history.csv/png · metrics_organ.csv · 그림 · pred.npz)
  fold##_<genus>/         : metrics_folds.csv - 내부 폴드별 1줄 + mean · median · std 줄
  ./                      : metrics_total.csv - genus별 1줄 + genus 간 mean · median · std 줄
                            model_final.pt · history_final.csv/png - 최종 모델
"""
import os, time
import numpy as np, pandas as pd, torch
from .batching import make_feeder
from .features import build_dataset
from .folds import nn_sim, nested_split
from .models import ds_blocks
from .plots import plot_history
from .train import train_model


def hms(sec):
    """초 -> 'H:MM:SS'"""
    s = int(sec)
    return f'{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}'


def fold_done(jdir, artifacts):
    """내부 폴드 산출물이 모두 있으면 True"""
    return all(os.path.exists(os.path.join(jdir, f)) and os.path.getsize(os.path.join(jdir, f)) > 0 for f in artifacts)


def with_stats(df, keys, stat_col):
    """표 df 에 수치 열의 mean · median · std 줄 3개를 붙인다 (keys = 그 줄에 채울 식별 열 dict, stat_col = 통계 이름 열)"""
    numc = df.select_dtypes(include=[np.number]).columns.drop([c for c in ('fold', 'inner') if c in df])
    agg  = pd.DataFrame([{**keys, stat_col: s, **{c: df[c].agg(s) for c in numc}} for s in ('mean', 'median', 'std')])
    return pd.concat([df, agg], ignore_index=True), agg


def run_cv(task, D, S, cfg, inputs, out_dir, name=''):
    """입력 조합 inputs 의 중첩 교차검증 -> genus별 결과 DataFrame [n_fold, *] (metrics_total.csv 도 저장)
       내부 폴드마다 데이터셋 생성(시험 · 검증 genus 를 참조에서 제외, 훈련 샘플 기준 표준화) -> 훈련(조기 종료) -> 테스트 genus 평가"""
    blocks = ds_blocks(inputs)
    feeder, lens = make_feeder(D, inputs, cfg.device, cfg.t_align)
    sc, W = task.score, max(len(g) for g in S.GENERA)
    tagw, pre = len(f'inner{cfg.inner_k:02d}'), f'{name} ' if name else ''
    rows_tot = []                                               # genus별 결과 [n_fold]

    for k in range(S.n_fold):                                   # 외부 폴드 k = 테스트 genus
        G, te = S.GENERA[k], S.GSMP[k]                          # 테스트 genus · 샘플 [n_test]
        fdir = os.path.join(out_dir, f'fold{k+1:02d}_{G}')
        os.makedirs(fdir, exist_ok=True)
        n_skip = sum(cfg.resume and fold_done(os.path.join(fdir, f'inner{j+1:02d}'), task.artifacts) for j in range(cfg.inner_k))
        print(f'\n=== {pre}fold {k+1}/{S.n_fold}  {G} === test {len(te):,} samples '
              f'({D.lab.org_id.iloc[te].nunique()} org{task.test_msg(te)})' + (f' | {n_skip} cached' if n_skip else ''))

        t_fold, rows_fold = time.perf_counter(), []
        for j in range(cfg.inner_k):                            # 내부 폴드 j = 검증 genus VAL_G[k][j]
            v, jtag = S.VAL_G[k][j], f'inner{j+1:02d}'
            jd = os.path.join(fdir, jtag)
            os.makedirs(jd, exist_ok=True)
            P = {a: os.path.join(jd, a) for a in task.artifacts}   # 산출물 경로
            idx = nested_split(S, k, j)

            cached = cfg.resume and fold_done(jd, task.artifacts)
            if cached:
                h = pd.read_csv(P['history.csv'])
                om, extra = task.load_cached(P)
                bep = int(h.epoch[h.is_best].iloc[0]) if h.is_best.any() else 0
            else:
                feeder.set_dataset(None)                        # 이전 데이터셋을 먼저 해제
                feeder.set_dataset(build_dataset(D, S, blocks, k, v, idx['train'], cfg.topk, cfg.t_eval, cfg.device))
                model, net, h, bep = train_model(task, inputs, feeder, lens, idx, cfg.seed + j,
                                                 f'{pre}{jtag:<{tagw}}', cfg, task.sample_weight(idx['train']))
                om, extra = task.evaluate(task.predict(net, feeder, lens, te), te, P,
                                          f'{pre}fold{k+1:02d} {G} {jtag} (val {S.GENERA[v]})')
                om.insert(0, 'fold', k + 1); om.insert(1, 'genus', G); om.insert(2, 'inner', j + 1)
                torch.save(model.state_dict(), P['model.pt'])
                h.round(cfg.csv_round).to_csv(P['history.csv'], index=False)
                plot_history(h, bep, P['history.png'], sc, task.loss_label, task.score_label)
                om.round(cfg.csv_round).to_csv(P['metrics_organ.csv'], index=False)
                del model, net
                if cfg.device != 'cpu':
                    torch.cuda.empty_cache()

            n_org, rest = task.fold_row(om, extra)
            r = dict(fold=k + 1, genus=G, inner=j + 1, val_genus=S.GENERA[v], val_sim=S.GSIM[k, v],
                     val_nn=nn_sim(S.GSIM, v, [k]), n_val=len(idx['val']), n_org=n_org, best_epoch=bep,
                     **{f'val_{sc}': float(h[f'vl_{sc}'].iloc[bep - 1]) if bep else np.nan}, **rest)
            rows_fold.append(r)
            print(f'  {jtag:<{tagw}} val {S.GENERA[v]:<{W}} (sim {S.GSIM[k, v]:+.2f}, nn {nn_sim(S.GSIM, v, [k]):.2f}) '
                  f'best@{bep:<4} val_{sc}={r[f"val_{sc}"]:.3f} {task.inner_msg(r)}' + (' | cached' if cached else ''))

        # metrics_folds.csv - 내부 폴드별 1줄 + mean · median · std 줄
        inn = pd.DataFrame(rows_fold)
        with_stats(inn, {'fold': k + 1, 'genus': G}, 'inner')[0].round(cfg.csv_round).to_csv(
            os.path.join(fdir, 'metrics_folds.csv'), index=False)

        # genus 성능 = 내부 폴드 지표의 평균 · 표준편차
        rows_tot.append(dict(fold=k + 1, genus=G, val_genera=';'.join(S.GENERA[S.VAL_G[k]]), n_test=len(te),
                             **task.genus_extra(te), n_org=int(inn.n_org.iloc[0]),
                             best_epoch_med=inn.best_epoch.median(), best_epoch_mean=inn.best_epoch.mean(),
                             best_epoch_sd=inn.best_epoch.std(),
                             **{kv: val for c in task.sd_cols for kv, val in ((c, inn[c].mean()), (c + '_sd', inn[c].std()))}))
        m  = rows_tot[-1]
        pm = lambda c: f'{m[c]:.3f}±{m[c + "_sd"]:.3f}'
        print(f'  {pre}fold {k+1}/{S.n_fold} {G} ({hms(time.perf_counter() - t_fold)}) '
              f'| best@{m["best_epoch_med"]:.1f}±{m["best_epoch_sd"]:.1f} val_{sc}={pm(f"val_{sc}")} {task.fold_msg(pm)}')

    # metrics_total.csv - genus별 1줄 + genus 간 mean · median · std 줄 (NaN genus 는 평균에서 빠진다)
    res = pd.DataFrame(rows_tot)
    with_stats(res, {'genus': ''}, 'fold')[0].round(cfg.csv_round).to_csv(os.path.join(out_dir, 'metrics_total.csv'), index=False)
    feeder.set_dataset(None)
    return res


def train_final(task, D, S, cfg, inputs, out_dir, res, name=''):
    """전체 샘플로 최종 모델 훈련 -> model (산출물이 이미 있고 resume 이면 None)
       epoch 수 = genus별 best epoch(내부 폴드 중앙값)들의 중앙값을 반올림
       문맥 특성은 자기 genus 만 참조에서 빼고, 표준화 기준은 전체 샘플이다"""
    n_ep = int(res.best_epoch_med.median() + 0.5)
    files = [os.path.join(out_dir, f) for f in ('model_final.pt', 'history_final.csv', 'history_final.png')]
    print(f'{name + " " if name else ""}final training epochs = {n_ep} (median of genus medians {res.best_epoch_med.median():.1f} | '
          f'mean of genus means {res.best_epoch_mean.mean():.1f})')
    if cfg.resume and all(os.path.exists(f) for f in files):
        print(f'  found {files[0]} -> skip')
        return None
    feeder, lens = make_feeder(D, inputs, cfg.device, cfg.t_align)
    allsmp = np.arange(task.n)                                  # 전체 샘플 번호 [n]
    feeder.set_dataset(build_dataset(D, S, ds_blocks(inputs), -1, -1, allsmp, cfg.topk, cfg.t_eval, cfg.device))
    model, net, hist, _ = train_model(task, inputs, feeder, lens, {'train': allsmp, 'val': allsmp[:0]},
                                      cfg.seed, f'{name} final'.strip(), cfg, task.sample_weight(allsmp), n_epoch=n_ep)
    torch.save(model.state_dict(), files[0])
    hist.drop(columns=['vl_loss', f'vl_{task.score}', 'is_best']).round(cfg.csv_round).to_csv(files[1], index=False)
    plot_history(hist, 0, files[2], task.score, task.loss_label, task.score_label)
    feeder.set_dataset(None)
    del net
    print(f'  saved {files[0]} | history_final.csv/png')
    return model


def print_summary(task, res, n_fold, inner_k, elapsed):
    """genus별 '내부 폴드 평균±표준편차' 표 + genus 간 mean · median · std + 성능 추정 + 상하위 genus"""
    _, agg = with_stats(res, {'genus': ''}, 'fold')
    W, MW, cols, tgt = max(map(len, res.genus)), 12, task.summary_cols, task.target
    head = f'{"#":>3}  {"genus":<{W}}  {"n":>7}  {"best@":^10}  ' + '  '.join(f'{c:^{MW}}' for c in cols)
    n_ok = int(res[tgt].notna().sum())                                      # 지표가 계산된 genus 수
    pmv  = lambda r, c: f'{r[c]:>6.3f}±{r[c + "_sd"]:<5.3f}'

    print(f'\nSUMMARY  {n_fold} genus x {inner_k} inner folds | total elapsed {hms(elapsed)}')
    print(f'{"=" * len(head)}\n{head}\n{"-" * len(head)}')
    for _, r in res.iterrows():
        print(f'{int(r.fold):>3}  {r.genus:<{W}}  {int(r.n_test):>7,}  {r.best_epoch_med:>5.1f}±{r.best_epoch_sd:<4.1f}  '
              + '  '.join(pmv(r, c) for c in cols))
    print('-' * len(head))
    for _, s in agg.iterrows():
        print(f'{"":>3}  {s.fold:<{W}}  {"":>7}  {s.best_epoch_med:>5.1f}{"":>5}  ' + '  '.join(f'{s[c]:>6.3f}{"":>6}' for c in cols))
    print('=' * len(head))
    print(f'± = sd across inner folds | mean/median/std = across {n_ok}/{n_fold} genus with computable metrics')

    print(f'\nPERFORMANCE ESTIMATE  1 vote per genus (n={n_ok})')
    print('  ' + ' | '.join(f'{c} {res[c].mean():.3f}±{res[c].std():.3f}' for c in task.perf_cols))
    print(f'  genus below {tgt} 0.1: {int((res[tgt] < 0.1).sum())}')
    print(f'\nper-genus ({tgt} descending, top/bottom 5)')
    top = res.sort_values(tgt, ascending=False)[['genus', 'n_test'] + task.top_cols]
    print((top if len(top) <= 10 else pd.concat([top.head(5), top.tail(5)])).round(3).to_string(index=False))


def summarize_ablations(task, RES, ablations, full_inputs, full_dir, out_dir, block_dim, genera, csv_round):
    """ablation별 genus 결과 RES {이름: DataFrame} -> ablation_summary.csv · genus x ablation 주 지표 표 저장 -> 요약 DataFrame
       full_dir/metrics_total.csv(전체 입력 모델)가 있으면 비교 기준 'full' 로 싣고 차이(delta_col) · 나빠진 genus 수를 더한다"""
    tgt, tabs = task.target, {}
    f_full = os.path.join(full_dir, 'metrics_total.csv')
    if os.path.exists(f_full):
        full = pd.read_csv(f_full)
        tabs['full'] = (full_inputs, full[pd.to_numeric(full.fold, errors='coerce').notna()])   # mean · median · std 줄 제외
    tabs.update({n: (ablations[n], r) for n, r in RES.items()})

    summ = pd.DataFrame([dict(ablation=n, inputs=' + '.join(inp), width=sum(block_dim[b] for b in inp),
                              n_genus=int(r[tgt].notna().sum()),
                              best_epoch_med=r.best_epoch_med.median(), best_epoch_mean=r.best_epoch_mean.mean(),
                              n_epoch_final=int(r.best_epoch_med.median() + 0.5),
                              **{kv: val for c in task.abl_cols for kv, val in ((c, r[c].mean()), (c + '_sd', r[c].std()))})
                         for n, (inp, r) in tabs.items()])
    gm = pd.DataFrame({n: r.set_index('genus')[tgt] for n, (_, r) in tabs.items()}).loc[genera]   # [n_fold, n_ablation]
    if 'full' in tabs:
        summ[task.delta_col] = summ[tgt] - summ[tgt][summ.ablation == 'full'].iloc[0]
        summ['n_genus_worse'] = [int((gm[n] < gm['full']).sum()) for n in summ.ablation]
    summ.round(csv_round).to_csv(os.path.join(out_dir, 'ablation_summary.csv'), index=False)
    gm.round(csv_round).to_csv(os.path.join(out_dir, task.genus_file), index_label='genus')

    cols = ['ablation', 'width', 'n_genus', 'n_epoch_final'] + task.abl_show \
           + ([task.delta_col, 'n_genus_worse'] if 'full' in tabs else [])
    print(summ[cols].round(3).to_string(index=False))
    print(f'\nsaved {out_dir}/ablation_summary.csv | {task.genus_file} | <ablation>/metrics_total.csv · model_final.pt')
    return summ
