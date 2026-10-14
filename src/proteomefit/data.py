"""base data(균주별 TSV) · 임베딩 h5 적재"""
import glob, os
from types import SimpleNamespace
import h5py, numpy as np, pandas as pd, torch
from tqdm import tqdm


def split_basedata(tsv, out_dir):
    """통합 TSV -> 균주별 TSV '순번_org_id.tsv' (행 · 문자열 그대로, 순번 = TSV 등장 순서), 반환 = 파일 수"""
    os.makedirs(out_dir, exist_ok=True)
    out, cur, i = None, None, 0
    with open(tsv, encoding='utf-8', newline='') as f:
        header = f.readline()
        for line in f:
            o = line.split('\t', 1)[0]                                  # 첫 열 = org_id
            if o != cur:                                                # 균주가 바뀌면 새 파일
                if out:
                    out.close()
                i, cur = i + 1, o
                out = open(os.path.join(out_dir, f'{i:02d}_{o}.tsv'), 'w', encoding='utf-8', newline='')
                out.write(header)
            out.write(line)
    out.close()
    return i


ID_DTYPE = {'org_id': str, 'org_name': str, 'gene_id': str, 'gene_name': str}   # ID 열은 문자열


def load_basedata(base_dir, **kw):
    """균주별 TSV 폴더 -> 통합 DataFrame [n_all, *] (파일명 순 = 원래 행 순서, kw 는 pd.read_csv 인자)"""
    kw.setdefault('dtype', ID_DTYPE)
    files = sorted(glob.glob(os.path.join(base_dir, '*.tsv')))
    return pd.concat([pd.read_csv(f, sep='\t', **kw) for f in files], ignore_index=True)


def load_samples(base_dir, label_col, finite_only):
    """라벨 · 메타 적재 -> 네임스페이스 D
       finite_only=False: fit_mean 공란(미상)만 제외 (분류, n 195,354)
       finite_only=True : 필수(-inf) · 미상 제외 (회귀, n 177,317)"""
    base = load_basedata(base_dir, usecols=['org_id', 'org_name', 'gene_id', 'seq_aa', 'seq_dna', label_col, 'comb_t'])
    base['genus'] = base.org_name.str.split().str[0]                    # org_name 첫 낱말 = genus
    yb   = base[label_col].to_numpy(np.float64)                         # 적합도 [n_all] (-inf = 필수, NaN = 미상)
    keep = np.flatnonzero(np.isfinite(yb) if finite_only else ~np.isnan(yb))   # 사용 샘플의 TSV(=h5) 번호 [n]
    lab  = base.iloc[keep].reset_index(drop=True)

    ocat = pd.Categorical(base.org_id)
    org_code_all = ocat.codes.astype(np.int32)                          # 균주 번호 [n_all]
    n_org = len(ocat.categories)                                        # 46

    # 균주 내 적합도 순위 [-1, 1] - 문맥 특성 otr_* 용 [n_all] (적합도 수치가 없으면 NaN)
    y_rank_all = np.full(len(base), np.nan, np.float32)
    for o in range(n_org):
        m = np.flatnonzero((org_code_all == o) & np.isfinite(yb))
        y_rank_all[m] = (np.argsort(np.argsort(yb[m])) + 0.5) / len(m) * 2 - 1

    return SimpleNamespace(
        base=base, lab=lab, keep=keep, yb=yb,
        ctb=base.comb_t.to_numpy(np.float64),                           # comb_t [n_all]
        org=lab.org_id.to_numpy(), gen=lab.genus.to_numpy(),            # 균주 · genus [n]
        org_code_all=org_code_all, org_code=org_code_all[keep], n_org=n_org,
        org2genus=pd.Series(base.genus.to_numpy()).groupby(org_code_all).first().to_numpy(),   # 균주 -> genus 이름 [n_org]
        y_rank_all=y_rank_all)


def ess_labels(D):
    """load_samples(finite_only=False) 결과 -> 필수 라벨 [n] (1 = 필수, fit_mean 이 -inf)"""
    return np.isneginf(D.yb[D.keep]).astype(np.int64)


def fit_labels(D):
    """load_samples(finite_only=True) 결과 -> 균주별 z-정규화 fit_mean [n]"""
    g = pd.Series(D.yb[D.keep]).groupby(D.org)
    return ((D.yb[D.keep] - g.transform('mean')) / g.transform('std')).to_numpy()


def load_embeddings(prot_h5, mean_h5, D):
    """ProteomeLM · ESM 토큰 평균 임베딩 -> (prot_raw [n_all, 1152], esm_mean_raw [n_all, 1152]) float32"""
    with h5py.File(prot_h5, 'r') as f:
        assert int(f.attrs['n_proteins']) == len(D.base), 'proteome h5/TSV size mismatch'
        assert np.array_equal(f['gene_id'][:].astype(str)[D.keep], D.lab.gene_id.to_numpy()), 'h5 samples misaligned'
        prot_raw = f['emb'][:].astype(np.float32)
    with h5py.File(mean_h5, 'r') as f:
        assert int(f.attrs['n_proteins']) == len(D.base), 'ESM mean h5/TSV size mismatch'
        esm_mean_raw = f['emb'][:].astype(np.float32)
    return prot_raw, esm_mean_raw


def load_esm_tokens(esm_h5, D, device, dtype, chunk, reserve_gib):
    """사용 샘플의 ESM 토큰을 GPU 한 배열에 이어 붙인다 -> (esm [총 토큰 수, 1152], esm_off [n + 1], esm_len [n])
       샘플 i = esm[esm_off[i]:esm_off[i+1]]. h5 는 chunk 개 단백질씩 읽는다"""
    n_all = len(D.base)
    with h5py.File(esm_h5, 'r') as f:
        assert int(f.attrs['n_proteins']) == n_all, 'ESM h5/TSV size mismatch'
        src, soff = f['emb'], f['offsets'][:].astype(np.int64)           # [총 토큰 수, 1152] fp16 · [n_all + 1]
        esm_len = (soff[D.keep + 1] - soff[D.keep]).astype(np.int64)     # 샘플별 토큰 수 [n]
        esm_off = np.concatenate([[0], np.cumsum(esm_len)])              # [n + 1]
        need, free = int(esm_off[-1]) * src.shape[1] * 2, torch.cuda.mem_get_info()[0]   # 필요 · 여유 VRAM (bytes)
        assert need + reserve_gib * 2**30 <= free, (f'ESM needs {need / 2**30:.1f} GiB + reserve {reserve_gib} GiB, '
                                                    f'free {free / 2**30:.1f} GiB')
        esm  = torch.empty((int(esm_off[-1]), src.shape[1]), dtype=dtype, device=device)
        used = np.zeros(n_all, bool)
        used[D.keep] = True                                              # TSV 행별 사용 여부 [n_all]
        dst = 0
        for a in tqdm(range(0, n_all, chunk), desc='ESM -> GPU'):
            z   = min(a + chunk, n_all)
            blk = src[int(soff[a]):int(soff[z])]                         # 이 묶음의 전 토큰 [t, 1152] fp16
            blk = blk[np.repeat(used[a:z], np.diff(soff[a:z + 1]))]      # 사용 샘플의 토큰만 [t', 1152]
            esm[dst:dst + len(blk)] = torch.from_numpy(blk).to(device)   # fp16 -> bf16
            dst += len(blk)
        assert dst == len(esm), 'ESM token count mismatch'
    torch.cuda.empty_cache()
    return esm, esm_off, esm_len
