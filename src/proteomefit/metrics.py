"""평가 지표 - 분류(필수 여부) · 회귀(적합도) · 신뢰도 보정, 균주별로 구해 평균한다"""
import numpy as np, pandas as pd, torch
from scipy.stats import rankdata
from sklearn.metrics import (confusion_matrix, matthews_corrcoef, f1_score, recall_score, precision_score,
                             accuracy_score, roc_auc_score, average_precision_score, log_loss)


def org_masks(org_ids, y_true=None):
    """균주 [m] -> 균주별 위치 목록. y_true [m] 를 주면 양·음성이 모두 있는 균주만"""
    out = []
    for o in pd.unique(org_ids):
        m = np.flatnonzero(org_ids == o)
        if y_true is None or 0 < y_true[m].sum() < len(m):
            out.append(m)
    return out


# ---------------------------------------------------------------- 분류

ESS_METRIC_COLS = ['tp', 'tn', 'fp', 'fn', 'mcc', 'f1', 'recall', 'precision', 'accuracy', 'roc_auc', 'pr_auc']
BOTH_CLASS_COLS = ['mcc', 'f1', 'recall', 'precision', 'roc_auc', 'pr_auc']   # 양·음성이 모두 필요한 지표


def mcc_org_mean(masks, y_true, logits, thresh):
    """균주별 MCC 의 평균 - 조기 종료 기준 (logits [m] GPU 텐서, masks = 그 안의 위치)"""
    pred = (torch.sigmoid(logits) >= thresh).int().cpu().numpy()     # [m] 0/1
    return float(np.mean([matthews_corrcoef(y_true[m], pred[m]) for m in masks]))


def ess_metrics(y_true, prob, thresh):
    """라벨 [m] · 양성 확률 [m] -> 지표 dict (키 = ESS_METRIC_COLS). 한 클래스만 있으면 BOTH_CLASS_COLS 는 NaN"""
    pred = (prob >= thresh).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred, labels=[0, 1]).ravel()
    out = dict(tp=int(tp), tn=int(tn), fp=int(fp), fn=int(fn),
               accuracy=accuracy_score(y_true, pred), **{c: np.nan for c in BOTH_CLASS_COLS})
    if 0 < y_true.sum() < len(y_true):
        out.update(mcc=matthews_corrcoef(y_true, pred),
                   f1=f1_score(y_true, pred, zero_division=0),
                   recall=recall_score(y_true, pred, zero_division=0),
                   precision=precision_score(y_true, pred, zero_division=0),
                   roc_auc=roc_auc_score(y_true, prob),
                   pr_auc=average_precision_score(y_true, prob))
    return {c: out[c] for c in ESS_METRIC_COLS}


def ess_metrics_by_organism(org_ids, y_true, prob, thresh):
    """균주 [m] · 라벨 [m] · 양성 확률 [m] -> 균주별 지표 DataFrame [n_org, 14] (org_id, n, n_pos + ESS_METRIC_COLS)"""
    out = []
    for o in pd.unique(org_ids):
        m = org_ids == o
        out.append(dict(org_id=o, n=int(m.sum()), n_pos=int(y_true[m].sum()), **ess_metrics(y_true[m], prob[m], thresh)))
    return pd.DataFrame(out).sort_values('org_id').reset_index(drop=True)


def in_total(om):
    """균주별 지표 DataFrame -> 총계에 넣을 균주 마스크 [n_org] (양·음성이 모두 있는 균주)"""
    return (om.n_pos > 0) & (om.n_pos < om.n)


# ---------------------------------------------------------------- 회귀

SIG_SUF  = '_sig'                                                   # 주 지표 접미사 (실측이 확실한 샘플, |comb_t| >= t_eval)
REG_COLS = ['spearman', 'pearson', 'ccc']
CLS_COLS = ['mcc', 'f1', 'recall', 'precision', 'accuracy', 'roc_auc', 'pr_auc', 'tp', 'tn', 'fp', 'fn']
FIT_METRIC_COLS = ['n' + SIG_SUF] + [m + s for s in (SIG_SUF, '_all') for m in REG_COLS + CLS_COLS]


def pearson(a, b):
    """a · b [m] -> Pearson"""
    return float(np.corrcoef(a, b)[0, 1])


def spearman(a, b):
    """Spearman = 순위의 Pearson"""
    return pearson(rankdata(a), rankdata(b))


def ccc(pred, obs):
    """예측 · 라벨 [m] -> 일치 상관계수(Concordance Correlation Coefficient, CCC)"""
    mp, mo = pred.mean(), obs.mean()
    return float(2 * np.mean((pred - mp) * (obs - mo)) / (pred.var() + obs.var() + (mp - mo) ** 2))


def sp_org_mean(masks, y_true, pred):
    """균주별 Spearman 의 평균 - 조기 종료 기준 (masks = y_true · pred [m] 안의 위치)"""
    return float(np.mean([spearman(pred[m], y_true[m]) for m in masks]))


def sign_metrics(y_true, mu):
    """라벨 · 예측 평균 [m] -> 0 기준 이진 분류 지표 dict (키 = CLS_COLS, 높음 = 1, AUC 점수 = mu)"""
    t, p = (y_true > 0).astype(int), (mu > 0).astype(int)          # 실제 · 예측 클래스 [m]
    tn, fp, fn, tp = confusion_matrix(t, p, labels=[0, 1]).ravel()
    both = 0 < t.sum() < len(t)                                     # 두 클래스가 모두 있는지
    return dict(mcc=matthews_corrcoef(t, p), f1=f1_score(t, p, zero_division=0),
                recall=recall_score(t, p, zero_division=0), precision=precision_score(t, p, zero_division=0),
                accuracy=accuracy_score(t, p), roc_auc=roc_auc_score(t, mu) if both else np.nan,
                pr_auc=average_precision_score(t, mu) if both else np.nan,
                tp=int(tp), tn=int(tn), fp=int(fp), fn=int(fn))


def fit_metrics(y_true, mu, sel):
    """라벨 · 예측 평균 · 실측이 확실한 샘플 마스크 (|comb_t| >= t_eval) [m] -> 지표 dict (키 = FIT_METRIC_COLS)
       _sig = 실측이 확실한 샘플 (주 지표), _all = 전 샘플"""
    out = {'n' + SIG_SUF: int(sel.sum())}
    for suf, a, b in ((SIG_SUF, mu[sel], y_true[sel]), ('_all', mu, y_true)):
        out.update({'spearman' + suf: spearman(a, b), 'pearson' + suf: pearson(a, b), 'ccc' + suf: ccc(a, b)})
        out.update({c + suf: val for c, val in sign_metrics(b, a).items()})
    return out


def fit_metrics_by_organism(org_ids, y_true, mu, sel):
    """균주 · 라벨 · 예측 평균 · 실측이 확실한 샘플 마스크 [m] -> 균주별 지표 DataFrame [n_org, 31] (org_id, n + FIT_METRIC_COLS)"""
    out = []
    for o in pd.unique(org_ids):
        m = org_ids == o
        out.append(dict(org_id=o, n=int(m.sum()), **fit_metrics(y_true[m], mu[m], sel[m])))
    return pd.DataFrame(out).sort_values('org_id').reset_index(drop=True)


# ---------------------------------------------------------------- 신뢰도 보정

PROB_COLS = ['ece', 'brier', 'log_loss', 'roc_auc']                       # 확률 지표
CALL_COLS = ['mcc', 'f1', 'recall', 'precision', 'accuracy', 'pr_auc']     # 확률 >= 0.5 판정 지표 (conf_ess)


def ece(hit, p, bins):
    """적중 [m] · 확률 [m] · 확률 구간 경계 -> 기대 보정 오차 (구간별 |확률 평균 - 적중 비율| 의 구간 크기 가중 평균)"""
    b = np.clip(np.digitize(p, bins) - 1, 0, len(bins) - 2)
    return float(sum((b == i).mean() * abs(p[b == i].mean() - hit[b == i].mean()) for i in np.unique(b)))


def prob_metrics(hit, p, bins):
    """적중 [m] · 확률 [m] -> 지표 dict (키 = PROB_COLS, 한 클래스만 있으면 roc_auc = NaN)"""
    return dict(ece=ece(hit, p, bins), brier=float(np.mean((p - hit) ** 2)), log_loss=log_loss(hit, p, labels=[0, 1]),
                roc_auc=roc_auc_score(hit, p) if 0 < hit.sum() < len(hit) else np.nan)


def calib_metrics_by_organism(org_ids, hit, probs, bins, call):
    """균주 · 적중 [m] · {단계: 확률 [m]} -> 균주별 지표 DataFrame [n_org, 3 + 지표 수 x 단계 수] (org_id, n, n_pos, {지표}_{단계})
       call = True 이면 확률 >= 0.5 판정 지표(CALL_COLS)도 구한다"""
    out = []
    for o in pd.unique(org_ids):
        m = org_ids == o
        r = dict(org_id=o, n=int(m.sum()), n_pos=int(hit[m].sum()))
        for stage, p in probs.items():
            met = prob_metrics(hit[m], p[m], bins)
            if call:
                e = ess_metrics(hit[m], p[m], 0.5)
                met.update({c: e[c] for c in CALL_COLS})
            r.update({f'{c}_{stage}': v for c, v in met.items()})
        out.append(r)
    return pd.DataFrame(out).sort_values('org_id').reset_index(drop=True)


def reliability_table(hit, probs, bins):
    """적중 [m] · {단계: 확률 [m]} -> 확률 구간별 표 [구간 수 x 단계 수, 5] (stage, bin, n, predicted, observed)"""
    out = []
    for stage, p in probs.items():
        t = pd.DataFrame({'p': p, 'hit': hit}).groupby(pd.cut(p, bins, include_lowest=True), observed=True).agg(
            n=('hit', 'size'), predicted=('p', 'mean'), observed=('hit', 'mean'))
        out.append(t.reset_index(names='bin').assign(stage=stage, bin=lambda d: d.bin.astype(str)))
    return pd.concat(out, ignore_index=True)[['stage', 'bin', 'n', 'predicted', 'observed']]
