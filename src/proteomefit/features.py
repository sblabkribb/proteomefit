"""입력 특성 - Biopython 서열 특성 · 정수 서열 · 직교상동체 목록 · 문맥 특성 · 데이터셋 (정의는 input_features.md)"""
from itertools import product
import numpy as np, pandas as pd, torch
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm import tqdm

AA     = list('ACDEFGHIKLMNPQRSTVWY')
AA_IDX = {a: i + 1 for i, a in enumerate(AA)}                 # 아미노산 -> 1 ~ 20
CODONS = [''.join(c) for c in product('TCAG', repeat=3)]      # 64코돈 고정 순서
CIDX   = {c: i for i, c in enumerate(CODONS)}
CTX_DIM   = 68                                                # 문맥 블록 폭
CTX_SPACE = {'ctx': 'esm', 'ctxp': 'prot', 'ctxb': 'both'}    # 문맥 블록 -> 유사도 공간


# ---------------------------------------------------------------- 서열 특성

def prot_feats(seq_aa):
    """단백질 서열 -> [72]: 물리화학 12 + 조성 20 + N말단 20잔기 조성 20 + C말단 20잔기 조성 20"""
    L = len(seq_aa)
    s = ''.join(ch for ch in seq_aa if ch in AA_IDX)          # ProtParam 은 표준 20종만 받는다
    if len(s) < 5:                                            # 너무 짧으면 길이만 채운다
        v = np.zeros(72, np.float32)
        v[0] = np.log1p(L)
        return v
    pa = ProteinAnalysis(s)
    h, t, e = pa.secondary_structure_fraction()
    cnt, n = pa.count_amino_acids(), len(s)
    comp = np.array([cnt[a] / n for a in AA], np.float32)                     # [20]
    nt = np.array([s[:20].count(a) / min(20, n) for a in AA], np.float32)     # [20]
    ct = np.array([s[-20:].count(a) / min(20, n) for a in AA], np.float32)    # [20]
    glob = np.array([
        np.log1p(L), pa.molecular_weight() / 1e4, pa.isoelectric_point(),
        pa.aromaticity(), pa.instability_index() / 100, pa.gravy(), h, t, e,
        pa.charge_at_pH(7.0) / 10,
        pa.molar_extinction_coefficient()[0] / 1e4,
        pa.molar_extinction_coefficient()[1] / 1e4], np.float32)             # [12]
    return np.concatenate([glob, comp, nt, ct])


def dna_feats(seq_dna):
    """CDS 서열 -> [68]: 코돈 빈도 64 + GC · GC1 · GC2 · GC3"""
    s = seq_dna.upper()
    n = len(s) // 3
    v = np.zeros(68, np.float32)
    if n < 3:
        return v
    cod = np.array([c for c in (CIDX.get(s[p:p + 3], -1) for p in range(0, n * 3, 3)) if c >= 0])
    if len(cod):
        v[:64] = np.bincount(cod, minlength=64) / len(cod)
    gc = np.isin(np.frombuffer(s.encode(), dtype='S1'), [b'G', b'C'])
    v[64] = gc.mean()
    for p in range(3):
        v[65 + p] = gc[p::3].mean() if len(gc[p::3]) else 0.0
    return v


def bio_features(seq_aa, seq_dna):
    """단백질 · CDS 서열 [n_all] -> Biopython 특성 [n_all, 140] float32"""
    return np.stack([np.concatenate([prot_feats(a), dna_feats(d)])
                     for a, d in tqdm(zip(seq_aa, seq_dna), total=len(seq_aa), desc='biopython')]).astype(np.float32)


def zscore_per_org(M, org_codes):
    """[n_all, d] -> 균주별 z-정규화 사본 [n_all, d] fp16"""
    Z = np.empty(M.shape, np.float16)
    for o in np.unique(org_codes):
        m = org_codes == o
        blk = M[m].astype(np.float32)
        Z[m] = ((blk - blk.mean(0)) / (blk.std(0) + 1e-6)).astype(np.float16)
    return Z


def encode_seqs(seqs, max_len):
    """단백질 서열 [n] -> 앞 max_len 잔기의 정수 인코딩을 한 배열에 이어 붙인 (flat [총 잔기 수] int8, off [n + 1], len [n], 원래 길이 [n])
       0 = 패딩, 1 ~ 20 = 아미노산, 21 = 그 외 (예: 'MKT' -> [11, 9, 17])"""
    lut = np.zeros(256, np.int8)
    for a, i in AA_IDX.items():
        lut[ord(a)] = i
    raw_len = np.fromiter(map(len, seqs), np.int64, len(seqs))
    seq_len = np.minimum(raw_len, max_len)
    seq_off = np.concatenate([[0], np.cumsum(seq_len)])
    flat    = np.empty(seq_off[-1], np.int8)
    for j, s in enumerate(seqs):
        u = lut[np.frombuffer(s[:max_len].encode(), np.uint8)]
        u[u == 0] = 21                                        # 표준 20종 밖의 문자
        flat[seq_off[j]:seq_off[j + 1]] = u
    return flat, seq_off, seq_len, raw_len


# ---------------------------------------------------------------- 직교상동체 목록

def whiten(E):
    """[n, d] -> 차원별 표준화 후 L2 정규화 [n, d] (내적 = 코사인)"""
    E = E - E.mean(0, keepdim=True)
    E = E / (E.std(0, keepdim=True) + 1e-6)
    return E / (E.norm(dim=1, keepdim=True) + 1e-6)


def org_slices(codes, n_org):
    """균주 번호 [n_all] -> 균주별 (start, stop) [n_org] (TSV 에서 균주가 연속 블록이어야 한다)"""
    bnd = np.flatnonzero(np.diff(codes) != 0) + 1
    span = np.concatenate([[0], bnd, [len(codes)]])
    sl = [None] * n_org
    for i in range(len(span) - 1):
        sl[codes[span[i]]] = (int(span[i]), int(span[i + 1]))
    assert all(s is not None for s in sl), 'organism blocks are not contiguous in the TSV'
    return sl


def build_ortho(E, slices, chunk, device, tag):
    """백색화 임베딩 [n_all, d] -> 균주별 최고 코사인 (sim [n_all, n_org] fp16, idx [n_all, n_org] int32)"""
    n, n_org = len(E), len(slices)
    Eg  = E.to(device)
    sim = np.empty((n, n_org), np.float16)
    idx = np.empty((n, n_org), np.int32)
    for s in tqdm(range(0, n, chunk), desc=f'ortho {tag}'):
        q = Eg[s:s + chunk]
        S = q @ Eg.T                                          # [c, n_all] 코사인
        S[torch.arange(len(q), device=device), torch.arange(s, s + len(q), device=device)] = -2.0   # 자기 자신 제외
        for o, (a, b) in enumerate(slices):
            v, j = S[:, a:b].max(1)
            sim[s:s + len(q), o] = v.cpu().numpy().astype(np.float16)
            idx[s:s + len(q), o] = (j + a).cpu().numpy()
        del S
    del Eg
    if device != 'cpu':
        torch.cuda.empty_cache()
    return sim, idx


def build_para(E, slices, k, device, tag):
    """백색화 임베딩 [n_all, d] -> 자기 게놈 안 상위 k 코사인 [n_all, k] fp16 (내림차순, 자기 자신 제외)"""
    out = np.empty((len(E), k), np.float16)
    for a, b in tqdm(slices, desc=f'para {tag}'):
        X = E[a:b].to(device)
        S = X @ X.T                                           # [m, m] 게놈 안 코사인
        S.fill_diagonal_(-2.0)
        out[a:b] = S.topk(k, 1).values.cpu().numpy().astype(np.float16)
        del S, X
    if device != 'cpu':
        torch.cuda.empty_cache()
    return out


def ortholog_lists(esm_mean_raw, prot_raw, org_code_all, n_org, k, chunk, device):
    """두 임베딩 [n_all, 1152] -> 유사도 공간 3종(esm · prot · both)의 직교상동체 목록 {공간: (sim, idx, para)}"""
    sl = org_slices(org_code_all, n_org)
    We, Wp = whiten(torch.from_numpy(esm_mean_raw)), whiten(torch.from_numpy(prot_raw))
    ort = {}
    for tag, E in (('esm', We), ('prot', Wp), ('both', torch.cat([We, Wp], 1) / np.sqrt(2.0))):
        ort[tag] = build_ortho(E, sl, chunk, device, tag) + (build_para(E, sl, k, device, tag),)
    return ort


# ---------------------------------------------------------------- 문맥 특성 · 데이터셋

def ctx_block(D, S, src, test_genus, val_genus, topk, t_eval):
    """유사도 공간 src · 시험 genus · 검증 genus -> (문맥 특성 [n, 68] float32, 이름 [68])
       참조 균주 = 시험 genus · 검증 genus · 질의 자신의 genus 어디에도 속하지 않는 균주"""
    sim_all, idx_all, para_all = D.ORT[src]
    keep, K = D.keep, max(topk)
    sim, idx = sim_all[keep].astype(np.float32), idx_all[keep]              # [n, 46]
    P  = para_all[keep].astype(np.float32)                                  # 자기 게놈 상위 코사인 [n, 20]
    n  = len(sim)
    og = pd.Categorical(D.org2genus, categories=S.GENERA).codes.astype(int) # 균주 -> genus 번호 [46]
    ok = (og[None, :] != S.gid[:, None]) & (og[None, :] != test_genus) & (og[None, :] != val_genus)   # 참조 가능 [n, 46]
    rbh_all = idx_all[idx, D.org_code[:, None]] == keep[:, None]            # 상호 최적 히트 [n, 46]

    # 참조 균주만 유사도 내림차순 상위 20개 -> 참조 단백질의 라벨 정보 [n, 20]
    smask = np.where(ok, sim, -np.inf)
    order = np.argsort(-smask, axis=1)[:, :K]
    r = np.arange(n)[:, None]
    S_, J  = smask[r, order], idx[r, order]                                 # 코사인 · 참조 단백질 번호
    valid  = np.isfinite(S_)
    yo     = D.yb[J].astype(np.float32)                                     # 참조 적합도
    ess    = np.isneginf(yo) & valid                                        # 참조가 필수
    fin    = np.isfinite(yo) & valid                                        # 참조 적합도 수치 있음
    strong = fin & (np.abs(np.nan_to_num(D.ctb[J])).astype(np.float32) >= t_eval)
    yo0    = np.where(fin, np.nan_to_num(yo), 0.0)
    yr     = D.y_rank_all[J]                                                # 참조 균주 내 순위
    finr   = np.isfinite(yr) & valid
    yr0    = np.where(finr, np.nan_to_num(yr), 0.0)
    Sf     = np.where(valid, S_, 0.0)
    R      = rbh_all[r, order] & valid

    feats, names = [], []
    def add(v, nm):
        feats.append(np.nan_to_num(np.asarray(v, np.float32), nan=0., posinf=0., neginf=0.))
        names.append(nm)
    def mean_over(mask, val, k):
        m = mask[:, :k]
        c = m.sum(1)
        return np.where(c > 0, np.where(m, val[:, :k], 0.).sum(1) / np.maximum(c, 1), 0.)

    with np.errstate(invalid='ignore', divide='ignore'):
        for k in topk:                                                  # cons_* 9열
            add(mean_over(valid, Sf, k), f'cons_simtop{k}')
        add(S_[:, 0] - S_[:, 2], 'cons_gap1_3')
        add(S_[:, 0] - S_[:, 9], 'cons_gap1_10')
        add(Sf.std(1), 'cons_sim_sd')
        add(valid.sum(1) / float(K), 'cons_nref')

        for k in topk:                                                  # otf_* 25열
            add(mean_over(fin, yo0, k), f'otf_mean{k}')
            add(ess[:, :k].sum(1) / np.maximum((fin | ess)[:, :k].sum(1), 1), f'otf_essfrac{k}')
            w = np.where(fin[:, :k], np.exp(5.0 * (Sf[:, :k] - 1.0)), 0.)
            add(np.where(w.sum(1) > 0, (w * yo0[:, :k]).sum(1) / np.maximum(w.sum(1), 1e-8), 0.),
                f'otf_wmean{k}')
        for k in topk:                                                  # otr_* 12열
            add(mean_over(finr, yr0, k), f'otr_mean{k}')
            w = np.where(finr[:, :k], np.exp(5.0 * (Sf[:, :k] - 1.0)), 0.)
            add(np.where(w.sum(1) > 0, (w * yr0[:, :k]).sum(1) / np.maximum(w.sum(1), 1e-8), 0.),
                f'otr_wmean{k}')
        add(np.where(finr[:, :10].any(1),
                     np.where(finr[:, :10], yr0[:, :10], np.inf).min(1), 0.), 'otr_min10')
        add((yr0[:, :20] < -0.5).sum(1) / 20., 'otr_frac_low20')
        add(mean_over(strong, yo0, 10), 'otf_mean_strong10')
        add(strong[:, :10].sum(1) / 10., 'otf_nstrong10')
        add(np.where(fin[:, :10].any(1),
                     np.where(fin[:, :10], yo0[:, :10], np.inf).min(1), 0.), 'otf_min10')
        add(np.where(fin[:, :10].any(1),
                     np.where(fin[:, :10], yo0[:, :10], -np.inf).max(1), 0.), 'otf_max10')
        add(np.where(fin[:, :10].any(1),
                     np.where(fin[:, :10], yo0[:, :10], 0.).std(1), 0.), 'otf_sd10')
        add(np.where(fin[:, 0], yo0[:, 0], 0.), 'otf_best_fit')
        add(S_[:, 0], 'otf_best_sim')
        add(ess[:, 0].astype(np.float32), 'otf_best_isess')
        add((yo0[:, :10] < -1).sum(1) / 10., 'otf_frac_lt_m1')
        add((yo0[:, :10] > 1).sum(1) / 10., 'otf_frac_gt_p1')

        add(R.sum(1) / float(K), 'rbh_frac')                            # rbh_* 14열
        add(mean_over(R, Sf, K), 'rbh_sim_mean')
        for k in (5, 10, 20):
            rf = R[:, :k] & fin[:, :k]
            add(rf.sum(1) / k, f'rbh_n{k}')
            add(mean_over(rf, yo0, k), f'rbh_fit{k}')
            add(mean_over(rf, yr0, k), f'rbh_rank{k}')
            add((R[:, :k] & ess[:, :k]).sum(1)
                / np.maximum((R[:, :k] & (fin | ess)[:, :k]).sum(1), 1), f'rbh_essfrac{k}')

        for k in topk:                                                  # para_* 8열 (라벨 미사용)
            add(P[:, :k].mean(1), f'para_simtop{k}')
        add(P[:, 0] - P[:, 2], 'para_gap1_3')
        add(P[:, 0] - P[:, 9], 'para_gap1_10')
        add(P.std(1), 'para_sim_sd')
    return np.stack(feats, 1).astype(np.float32), names


def build_dataset(D, S, blocks, k, v, ref, topk, t_eval, device):
    """데이터셋 블록 목록 · 시험 genus k · 검증 genus v · 표준화 기준 샘플 ref -> 데이터셋 [n, d] fp16 GPU 텐서
       예: blocks ['prot', 'bio', 'ctx', 'ctxp', 'ctxb'] -> [n, 1496]. k = v = -1 이면 자기 genus 만 참조에서 뺀다 (최종 모델)"""
    fixed = {'esm_mean': getattr(D, 'EMEAN', None), 'prot': D.PROT, 'bio': D.BIO}
    M = np.concatenate([fixed[b].astype(np.float32) if b in fixed
                        else ctx_block(D, S, CTX_SPACE[b], k, v, topk, t_eval)[0] for b in blocks], 1)
    M -= M[ref].mean(0)                                                     # 기준 샘플의 열별 평균 · 표준편차로 표준화
    M /= M[ref].std(0) + 1e-6
    return torch.from_numpy(np.nan_to_num(M, copy=False).astype(np.float16)).to(device)


def block_dims(D, seq_dim, esm_dim=1152):
    """입력별 연결 벡터 폭 dict (예: {'esm': 1152, 'seq': 64, 'prot': 1152, 'bio': 140, 'ctx': 68, ...})"""
    return {'esm': esm_dim, 'seq': seq_dim, 'esm_mean': esm_dim, 'prot': D.PROT.shape[1],
            'bio': D.BIO.shape[1], 'ctx': CTX_DIM, 'ctxp': CTX_DIM, 'ctxb': CTX_DIM}
