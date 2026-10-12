"""genus 단위 중첩 교차검증 분할 - 외부 폴드 = 시험 genus, 내부 폴드 = 외삽 거리가 비슷한 검증 genus"""
from types import SimpleNamespace
import numpy as np, pandas as pd
from sklearn.metrics.pairwise import cosine_similarity


def genus_order(genus):
    """genus 이름 [n] -> 샘플 수 내림차순 genus 이름 [n_genus] (외부 폴드 순서)"""
    return pd.Series(genus).value_counts().index.to_numpy()


def nn_sim(sim, g, excl):
    """genus 유사도 sim [G, G] · genus g · 제외 genus 목록 excl -> g 와 나머지 genus 중 최근접의 유사도 (외삽 거리, 클수록 가까운 친척이 있다)"""
    pool = [h for h in range(len(sim)) if h != g and h not in excl]
    return float(sim[g, pool].max())


def pick_val(k, sim, ok, n_val):
    """외부 폴드 k · genus 유사도 sim [G, G] · 후보 가능 여부 ok [G] · 개수 n_val -> 검증 genus 번호 [n_val]
       외삽 거리 차이 |nn_sim(v, {k}) - nn_sim(k, {})| 가 작은 순 (예: 시험 0.39 -> 후보 0.41 차이 0.02 가 0.89 차이 0.50 보다 먼저)"""
    rest = [g for g in range(len(sim)) if g != k]
    nn_k = rest[int(np.argmax(sim[k, rest]))]                      # 시험 genus 의 최근접 genus - 훈련 · 문맥 참조에 남긴다
    s_k  = nn_sim(sim, k, [])                                      # 시험 genus 의 외삽 거리
    cand = np.array([v for v in rest if v != nn_k and ok[v]])      # 검증 후보 [m]
    gap  = np.abs([nn_sim(sim, v, [k]) - s_k for v in cand])       # 후보별 외삽 거리 차이 [m]
    return cand[np.argsort(gap, kind='stable')[:n_val]]


def both_classes(y, org):
    """샘플 묶음의 라벨 y [m] · 균주 org [m] -> 양·음성이 모두 있는 균주가 하나 이상인지 (분류 검증 genus 후보 조건)"""
    return any(0 < y[org == o].sum() < (org == o).sum() for o in np.unique(org))


def make_splits(gen, org2genus, org_vec, inner_k, candidate=None):
    """샘플 genus [n] · 균주 -> genus [n_org] · 균주 프로테옴 표현 [n_org, d] -> 분할 네임스페이스 S
       candidate(genus 의 샘플 번호) -> 검증 genus 후보 여부 (None 이면 전부 후보)"""
    GENERA = genus_order(gen)                                           # 폴드 순서의 genus 이름 [n_fold]
    n_fold = len(GENERA)                                                # 28
    gid    = pd.Categorical(gen, categories=GENERA).codes.astype(int)   # 샘플 -> genus 번호 [n]
    GSMP   = [np.flatnonzero(gid == g) for g in range(n_fold)]          # genus별 샘플 번호 [n_fold][*]

    # genus 표현 -> 중심화 코사인 유사도 -> 검증 genus 선택
    org_gid = pd.Categorical(org2genus, categories=GENERA).codes        # 균주 -> genus 번호 [n_org]
    GVEC    = pd.DataFrame(org_vec).groupby(org_gid).mean().to_numpy()  # genus 표현 [n_fold, d]
    GSIM    = cosine_similarity(GVEC - GVEC.mean(0))                    # [n_fold, n_fold], 값 범위 [-1, 1]
    VAL_OK  = np.array([candidate(GSMP[g]) if candidate else True for g in range(n_fold)])   # 검증 genus 후보 여부 [n_fold]
    VAL_G   = [pick_val(k, GSIM, VAL_OK, inner_k) for k in range(n_fold)]                   # 검증 genus 번호 [n_fold][inner_k]
    return SimpleNamespace(GENERA=GENERA, n_fold=n_fold, gid=gid, GSMP=GSMP, GSIM=GSIM, VAL_OK=VAL_OK, VAL_G=VAL_G)


def nested_split(S, k, j):
    """외부 폴드 k · 내부 폴드 j -> {'train', 'val', 'test': 샘플 번호 (오름차순)}"""
    v = S.VAL_G[k][j]
    return {'train': np.flatnonzero((S.gid != k) & (S.gid != v)), 'val': S.GSMP[v], 'test': S.GSMP[k]}


def print_splits(S, lab, inner_k, y=None):
    """분할 요약 출력 + 누출 점검 (출력이 warn 없이 끝나면 정상). y 를 주면 양성 비율도 출력"""
    W = max(len(g) for g in S.GENERA)
    print(f'nested CV: {S.n_fold} outer x {inner_k} inner (genus-level) = {S.n_fold * inner_k:,} trainings / {len(lab):,} samples')
    print(f'val-genus candidates {int(S.VAL_OK.sum())}/{S.n_fold} (excluded: {", ".join(S.GENERA[~S.VAL_OK]) or "none"})\n')
    for k in range(S.n_fold):
        te  = S.GSMP[k]
        pos = f', pos {y[te].mean() * 100:4.1f}%' if y is not None else ''
        print(f'fold{k+1:02d} {S.GENERA[k]:<{W}}: test {len(te):>6,}({lab.org_id.iloc[te].nunique()}org{pos}) | val '
              + ' '.join(f'{S.GENERA[v]}({S.GSIM[k, v]:+.2f}, nn {nn_sim(S.GSIM, v, [k]):.2f}, {len(S.GSMP[v]):,})' for v in S.VAL_G[k]))

    for k in range(S.n_fold):
        if y is not None and y[S.GSMP[k]].sum() == 0:                   # 지표가 NaN 이 될 테스트 genus
            print(f'note fold{k+1:02d} {S.GENERA[k]}: no positive samples')
        for j in range(inner_k):
            d = nested_split(S, k, j)
            if (np.isin(S.gid[d['train']], [k, S.VAL_G[k][j]]).any() or (S.gid[d['val']] != S.VAL_G[k][j]).any()
                    or len(np.unique(S.VAL_G[k])) < inner_k):
                print(f'warn fold{k+1:02d}/inner{j+1:02d}: genus leak or duplicate val genus')
