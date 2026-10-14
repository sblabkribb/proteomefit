"""과제 정의 - 라벨 · 모델 · 손실 · 조기 종료 지표 · 평가 산출물 · 결과 표 형식

EssTask : 필수 유전자 이진 분류 (가중 BCE, 균주 평균 MCC 로 조기 종료)
FitTask : 적합도 회귀 (미니 MLP 1000개, 가지별 가중 1 - CCC, 균주 평균 Spearman 으로 조기 종료, 주 지표 = 실측이 확실한 샘플)
두 클래스는 같은 메서드를 가진다 (train.train_model · cv.run_cv 가 호출).
"""
import numpy as np, pandas as pd, torch, torch.nn as nn
from .batching import predict
from .data import ess_labels, fit_labels
from .folds import both_classes
from .metrics import (org_masks, mcc_org_mean, ess_metrics_by_organism, in_total, org_summary, ESS_METRIC_COLS,
                      ccc, sp_org_mean, fit_metrics_by_organism, FIT_METRIC_COLS, SIG_SUF)
from .models import EssNet, FitNet
from .plots import scatter_plot, sign_confusion_plot
from .train import taxon_weight, ccc_loss, bootstrap_counts


class EssTask:
    """필수 유전자 이진 분류 - 라벨 1 = 필수(fit_mean -inf), 0 = 비필수(수치)
       D: data.load_samples(finite_only=False) 결과 (D.block_dim 은 모델 생성 전에 채운다), cfg: config.Config"""
    score, target = 'mcc', 'mcc'                                     # 조기 종료 지표 이름 · 결과 정렬 기준 열
    classes   = ('non-essential', 'essential')                       # 라벨 0 / 1
    artifacts = ('model.pt', 'history.csv', 'history.png', 'metrics_organ.csv', 'pred.npz')
    loss_label, score_label = 'weighted BCE loss', 'MCC (per-organism mean)'
    sd_cols      = ['val_mcc'] + ESS_METRIC_COLS                     # genus 표에 평균 · _sd 를 남길 열
    summary_cols = ['val_mcc', 'mcc', 'recall', 'precision', 'roc_auc', 'pr_auc']   # 요약 표 열
    top_cols     = summary_cols[1:]                                  # genus 상하위 표 열
    perf_cols    = ['mcc', 'f1', 'recall', 'precision', 'accuracy', 'roc_auc', 'pr_auc']   # 성능 추정 출력 열
    brief_cols   = ['mcc', 'roc_auc', 'pr_auc']                      # ablation 완료 한 줄 출력 열
    abl_cols     = perf_cols                                         # ablation_summary.csv 에 genus 평균 · _sd 를 남길 열
    abl_show     = ['mcc', 'mcc_sd', 'roc_auc', 'pr_auc']            # ablation 요약 출력 열
    delta_col, genus_file = 'd_mcc', 'ablation_genus_mcc.csv'

    def __init__(self, D, cfg):
        self.D, self.cfg = D, cfg
        self.y  = ess_labels(D)                                      # 라벨 [n] (1 = 필수)
        self.yg = torch.from_numpy(self.y).float().to(cfg.device)   # GPU 라벨 [n]
        self.n  = len(self.y)

    def val_candidate(self, smp):
        """검증 genus 후보 - 검증 샘플에 양·음성이 모두 있는 균주가 하나 이상"""
        return both_classes(self.y[smp], self.D.org[smp])

    def build_model(self, inputs, seed):
        c = self.cfg
        return EssNet(inputs, self.D.block_dim, mlp_hidden=c.mlp_hidden, droprate=c.droprate,
                      seq_dim=c.seq_dim, seq_head=c.seq_head, n_head=c.n_head, max_len=c.max_len)

    def arch_msg(self, model):
        return f'MLP {" -> ".join(map(str, [model.d_in, *self.cfg.mlp_hidden, 1]))}'

    def sample_weight(self, smp):
        """훈련 샘플 번호 [m] -> 샘플 가중치 [m] (평균 1) = 분류군 역크기"""
        w = taxon_weight(self.D.gen[smp], self.D.org[smp])
        return (w / w.mean()).astype(np.float32)

    # ---- 훈련 (train.train_model 이 호출)
    def begin(self, tr, sw_full, seed):
        """손실 준비 - 양성 가중 = 훈련 샘플의 음성 수 / 양성 수"""
        p = self.y[tr].sum()
        pw = torch.tensor((len(tr) - p) / max(p, 1), dtype=torch.float32, device=self.cfg.device)
        self.lossf, self.sw = nn.BCEWithLogitsLoss(pos_weight=pw, reduction='none'), sw_full

    def batch_loss(self, out, b):
        """로짓 [B] · GPU 샘플 번호 b [B] -> 가중 BCE (스칼라)"""
        w = self.sw[b]
        return (self.lossf(out, self.yg[b]) * w).sum() / w.sum()

    def point(self, out):
        return out                                                  # 로짓 [B]

    def val_samples(self, vl):
        return vl                                                   # 검증 genus 전 샘플

    def predict(self, net, feeder, lens, smp):
        """-> smp 순서의 로짓 [m] fp32 GPU 텐서"""
        c = self.cfg
        return predict(net, feeder, lens, smp, c.eval_batch, c.eval_tok_budget, c.amp_dtype)

    def val_point(self, raw):
        return raw

    def val_loss(self, raw, smp):
        return self.lossf(raw, self.yg[torch.from_numpy(smp).to(self.cfg.device)]).mean().item()   # 클래스 가중만

    def scorer(self, smp):
        """샘플 번호 [m] -> (smp 순서 로짓 [m] -> 균주 평균 MCC) 함수 (양·음성이 모두 있는 균주만)"""
        masks, ys = org_masks(self.D.org[smp], self.y[smp]), self.y[smp]
        return lambda p: mcc_org_mean(masks, ys, p, self.cfg.thresh)

    # ---- 평가 · 결과 표 (cv.run_cv 가 호출)
    def evaluate(self, raw, te, P, title):
        """테스트 로짓 [n_test] -> (균주별 지표 DataFrame, 추가 값 dict). pred.npz 저장 (혼동행렬은 calib.nested_cv)"""
        pb = torch.sigmoid(raw).cpu().numpy()                       # 양성 확률 [n_test]
        yt, th = self.y[te], self.cfg.thresh
        om = ess_metrics_by_organism(self.D.org[te], yt, pb, th)
        np.savez(P['pred.npz'], prob=pb)                            # [n_test] - 테스트 genus 샘플 순
        return om, {}

    def load_cached(self, P):
        return pd.read_csv(P['metrics_organ.csv']), {}

    def fold_row(self, om, extra):
        """균주별 지표 -> (지표가 계산된 균주 수, 균주 평균 · 표준편차 dict) - 양·음성이 모두 있는 균주만"""
        tot = om[in_total(om)]
        return len(tot), org_summary(tot, ESS_METRIC_COLS)

    def genus_extra(self, te):
        return {'n_pos': int(self.y[te].sum())}

    def test_msg(self, te):
        return f', pos {self.y[te].mean() * 100:.1f}%'

    def inner_msg(self, r):
        return (f'| test mcc={r["mcc"]:.3f} recall={r["recall"]:.3f} precision={r["precision"]:.3f} '
                f'roc_auc={r["roc_auc"]:.3f} pr_auc={r["pr_auc"]:.3f}')

    def fold_msg(self, pm):
        return (f'| mcc={pm("mcc")} recall={pm("recall")} precision={pm("precision")} '
                f'roc_auc={pm("roc_auc")} pr_auc={pm("pr_auc")}')


class FitTask:
    """적합도 회귀 - 라벨 = 균주별 z-정규화 fit_mean, 예측 = 미니 MLP 가지 예측의 평균 mu · 표준편차 sd
       D: data.load_samples(finite_only=True) 결과, cfg: config.Config"""
    score, target = 'sp', 'spearman' + SIG_SUF
    classes   = ('low', 'high')                                      # 0 기준 이진 분류 (라벨 0 / 1)
    artifacts = ('model.pt', 'history.csv', 'history.png', 'metrics_organ.csv', 'scatter.png', 'confusion.png', 'pred.npz')
    loss_label, score_label = '1 - CCC loss', 'Spearman (|comb_t| >= 2)'
    show         = [m + s for s in (SIG_SUF, '_all') for m in ('spearman', 'ccc', 'mcc', 'roc_auc')]
    sd_cols      = ['val_sp'] + FIT_METRIC_COLS
    summary_cols = ['val_sp'] + show
    top_cols = perf_cols = brief_cols = show
    abl_cols     = FIT_METRIC_COLS
    abl_show     = show + ['spearman' + SIG_SUF + '_sd']
    delta_col, genus_file = 'd_sp', 'ablation_genus_spearman.csv'
    val_candidate = None                                             # 검증 genus 후보 = 시험 genus 와 그 최근접을 뺀 전부

    def __init__(self, D, cfg):
        self.D, self.cfg = D, cfg
        self.y   = fit_labels(D)                                     # 라벨 = 균주별 z-정규화 fit_mean [n]
        self.tt  = np.abs(np.nan_to_num(D.ctb[D.keep]))              # |comb_t| [n] (훈련 샘플 가중치)
        self.ext = self.tt >= cfg.t_eval                             # |comb_t| >= 2 [n] (조기 종료 검증 샘플)
        self.yg  = torch.from_numpy(self.y).float().to(cfg.device)
        self.n   = len(self.y)

    def build_model(self, inputs, seed):
        c = self.cfg
        return FitNet(inputs, self.D.block_dim, seed, n_branch=c.n_branch, br_feat=c.br_feat, br_hidden=c.br_hidden,
                      block_p=c.block_p, droprate=c.droprate, seq_dim=c.seq_dim, seq_head=c.seq_head,
                      n_head=c.n_head, max_len=c.max_len)

    def arch_msg(self, model):
        c = self.cfg
        return f'concat {model.d_in} -> {c.n_branch} x MLP({" -> ".join(map(str, [c.br_feat, *c.br_hidden, 1]))}), block p={c.block_p}'

    def sample_weight(self, smp):
        """훈련 샘플 번호 [m] -> 샘플 가중치 [m] (평균 1) = clip(|comb_t|) x 분류군 역크기"""
        w = taxon_weight(self.D.gen[smp], self.D.org[smp]) * np.clip(self.tt[smp], *self.cfg.w_clip)
        return (w / w.mean()).astype(np.float32)

    # ---- 훈련
    def begin(self, tr, sw_full, seed):
        """손실 준비 - 가지별 부트스트랩 등장 횟수 [n, n_branch] (bootstrap=True 일 때)"""
        c = self.cfg
        self.sw  = sw_full
        self.cnt = bootstrap_counts(tr, self.n, c.n_branch, seed, c.device) if c.bootstrap else None

    def batch_loss(self, out, b):
        """가지별 예측 [B, K] · GPU 샘플 번호 b [B] -> 가지별 가중 1 - CCC 의 평균"""
        w = self.sw[b][:, None] * self.cnt[b] if self.cnt is not None else self.sw[b]   # [B, K] 또는 [B]
        return ccc_loss(out, self.yg[b], w)

    def point(self, out):
        return out.mean(1)                                          # 예측 평균 [B]

    def val_samples(self, vl):
        return vl[self.ext[vl]]                                     # 검증 genus 중 |comb_t| >= 2

    def predict(self, net, feeder, lens, smp):
        """-> smp 순서의 (예측 평균 mu [m], 가지 간 표준편차 sd [m]) numpy"""
        c = self.cfg
        r = predict(net, feeder, lens, smp, c.eval_batch, c.eval_tok_budget, c.amp_dtype,
                    reduce=lambda P: torch.stack([P.mean(1), P.std(1)], 1))   # [m, 2]
        return r[:, 0].cpu().numpy(), r[:, 1].cpu().numpy()

    def val_point(self, raw):
        return raw[0]

    def val_loss(self, raw, smp):
        return 1 - ccc(raw[0], self.y[smp])

    def scorer(self, smp):
        """샘플 번호 [m] -> (smp 순서 예측 평균 [m] -> |comb_t| >= 2 샘플의 균주 평균 Spearman) 함수"""
        xt = self.ext[smp]
        masks, ys = org_masks(self.D.org[smp][xt]), self.y[smp][xt]
        def f(p):
            p = p.cpu().numpy() if torch.is_tensor(p) else p
            return sp_org_mean(masks, ys, p[xt])
        return f

    # ---- 평가 · 결과 표
    def evaluate(self, raw, te, P, title):
        """테스트 (mu, sd) -> (균주별 지표 DataFrame, {}). scatter.png · confusion.png · pred.npz 저장
           주 지표 _sig = 실측이 확실한 샘플 (|comb_t| >= t_eval), 보조 지표 _all = 전 샘플"""
        mu, sd = raw
        yt, sel = self.y[te], self.ext[te]                          # 라벨 · 실측이 확실한 샘플 마스크 [n_test]
        om  = fit_metrics_by_organism(self.D.org[te], yt, mu, sel)
        nm_sig = f'|comb_t| >= {self.cfg.t_eval:g} (n={int(sel.sum()):,})'
        nm_all = f'all (n={len(te):,})'
        scatter_plot([(f'{nm_sig}\nSpearman {om["spearman" + SIG_SUF].mean():.3f}  CCC {om["ccc" + SIG_SUF].mean():.3f}', yt[sel], mu[sel], 1.96 * sd[sel]),
                      (f'{nm_all}\nSpearman {om.spearman_all.mean():.3f}  CCC {om.ccc_all.mean():.3f}', yt, mu, 1.96 * sd)],
                     P['scatter.png'], title)
        sign_confusion_plot([(f'{nm_sig}\nMCC {om["mcc" + SIG_SUF].mean():.3f}', yt[sel], mu[sel]),
                             (f'{nm_all}\nMCC {om.mcc_all.mean():.3f}', yt, mu)], P['confusion.png'], title, self.classes)
        np.savez(P['pred.npz'], mu=mu, sd=sd)
        return om, {}

    def load_cached(self, P):
        return pd.read_csv(P['metrics_organ.csv']), {}

    def fold_row(self, om, extra):
        """균주별 지표 -> (주 지표가 계산된 균주 수, 균주 평균 · 표준편차 dict, n_sig 은 균주 합)"""
        rest = org_summary(om, FIT_METRIC_COLS)
        rest['n' + SIG_SUF] = int(om['n' + SIG_SUF].sum())
        return int(om[self.target].notna().sum()), rest

    def genus_extra(self, te):
        return {}

    def test_msg(self, te):
        return ''

    def inner_msg(self, r):
        return (f'| Sp_sig={r["spearman_sig"]:.3f} MCC_sig={r["mcc_sig"]:.3f} '
                f'| Sp_all={r["spearman_all"]:.3f} MCC_all={r["mcc_all"]:.3f}')

    def fold_msg(self, pm):
        return (f'| Sp_sig={pm("spearman_sig")} MCC_sig={pm("mcc_sig")} '
                f'| Sp_all={pm("spearman_all")} MCC_all={pm("mcc_all")}')

    def predict_table(self, net, feeder, lens, smp):
        """데이터셋이 준비된 샘플 번호 smp [m] -> DataFrame [m, 4] (gene_id, org_id, mu, sd)"""
        mu, sd = self.predict(net, feeder, lens, smp)
        return pd.DataFrame({'gene_id': self.D.lab.gene_id.to_numpy()[smp], 'org_id': self.D.org[smp], 'mu': mu, 'sd': sd})
