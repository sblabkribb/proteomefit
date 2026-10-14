"""배치 구성 · GPU 상주 배치 공급 · 추론"""
import numpy as np, torch


def make_batches(smp, lens, max_size, tok_budget, rng=None, pool_mult=50):
    """샘플 번호 [m] -> 배치 목록 (배치 = 샘플 번호 배열)
       |배치| <= max_size, |배치| x 배치 내 최대 길이 <= tok_budget
       rng 지정(훈련): 섞은 뒤 pool_mult*max_size 구간 안에서 길이순 정렬, 배치 순서도 섞는다
       rng=None(평가): 전체 길이순 정렬"""
    if rng is None:
        order = smp[np.argsort(lens[smp], kind='stable')]
    else:
        r = rng.permutation(smp)
        order = np.concatenate([c[np.argsort(lens[c], kind='stable')]
                                for c in np.array_split(r, max(len(r) // (pool_mult * max_size), 1))])
    out, cur, mx = [], [], 0
    for i in order:
        m = max(mx, int(lens[i]))
        if cur and (len(cur) + 1 > max_size or (len(cur) + 1) * m > tok_budget):
            out.append(np.array(cur)); cur = []
            m = int(lens[i])
        cur.append(i); mx = m
    out.append(np.array(cur))
    return [out[k] for k in rng.permutation(len(out))] if rng is not None else out


def gather_pad(flat, off, ln, b, lmax, t_align):
    """이어 붙인 GPU 배열 flat [총 길이, ...] · 시작 위치 off [n] · 길이 ln [n] · GPU 샘플 번호 b [B] · 배치 최대 길이 lmax
       -> (배치 [B, L, ...], 패딩 마스크 [B, L] True = 패딩), L = lmax 를 t_align 배수로 올린 값
       패딩 자리에는 이웃 샘플의 원소가 들어가며 마스크로 가린다"""
    L   = -(-lmax // t_align) * t_align
    pos = torch.arange(L, device=flat.device)                              # [L]
    at  = (off[b][:, None] + pos).clamp_(max=len(flat) - 1)                # [B, L] 원소 위치
    return flat[at], pos >= ln[b][:, None]


class BatchFeeder:
    """배치 공급기 - ESM 토큰 · 정수 서열 · 데이터셋이 모두 GPU 상주 텐서이며, 샘플 번호로 색인해 배치를 만든다
       esm · sflat 이 None 이면 그 입력은 만들지 않는다 (None 으로 넘긴다)"""

    def __init__(self, esm, eoff, elen, sflat, soff, slen, device, t_align=8):
        self.dev, self.ds, self.esm, self.sflat, self.t_align = device, None, esm, None, t_align   # esm = [총 토큰 수, 1152] bf16 GPU
        self.elen, self.slen = elen, slen                                  # 배치 최대 길이 계산용 numpy [n]
        if esm is not None:
            self.eoff, self.elen_g = (torch.from_numpy(a).to(device) for a in (eoff[:-1], elen))          # [n] · [n]
        if sflat is not None:
            self.sflat, self.soff, self.slen_g = (torch.from_numpy(a).to(device) for a in (sflat, soff[:-1], slen))   # int8 [총 잔기 수] · [n] · [n]

    def set_dataset(self, ds):
        """데이터셋 [n, d] fp16 GPU 텐서를 교체한다 (None 이면 해제)"""
        self.ds = ds

    def stream(self, batches):
        """배치 목록 -> (tok [B,T,1152] | None, pad [B,T] | None, seq [B,L] | None, ds [B,d], b [B] GPU 샘플 번호) 생성기"""
        idx = torch.from_numpy(np.concatenate(batches)).to(self.dev).split([len(bn) for bn in batches])   # 배치별 GPU 샘플 번호
        for bn, b in zip(batches, idx):
            tok = pad = seq = None
            if self.esm is not None:
                tok, pad = gather_pad(self.esm, self.eoff, self.elen_g, b, int(self.elen[bn].max()), self.t_align)
            if self.sflat is not None:
                seq, spad = gather_pad(self.sflat, self.soff, self.slen_g, b, int(self.slen[bn].max()), self.t_align)
                seq = seq.masked_fill(spad, 0).long()                     # 서열 패딩 = 0
            yield tok, pad, seq, self.ds[b], b


def make_feeder(D, inputs, device, t_align):
    """입력 조합 inputs 의 배치 공급기 · 배치 구성용 길이 [n] (ESM 토큰 수 > 서열 길이 > 1 우선)"""
    use_esm, use_seq = 'esm' in inputs, 'seq' in inputs
    lens = D.esm_len if use_esm else D.seq_len if use_seq else np.ones(len(D.lab), np.int64)
    feeder = BatchFeeder(D.esm if use_esm else None, D.esm_off, D.esm_len,
                         D.seq_flat if use_seq else None, D.seq_off, D.seq_len, device, t_align)
    return feeder, lens


def predict(net, feeder, lens, smp, max_size, tok_budget, amp_dtype, reduce=None):
    """eval 모드 추론 -> smp 순서의 출력 fp32 GPU 텐서 (분류 로짓 [m] / 회귀는 reduce 로 배치마다 줄인 값)"""
    bs = make_batches(smp, lens, max_size, tok_budget)
    net.eval(); out = []
    with torch.no_grad(), torch.autocast('cuda', dtype=amp_dtype):
        for *xb, ds, _ in feeder.stream(bs):
            o = net(*xb, ds).float()
            out.append(reduce(o) if reduce else o)
    back = torch.from_numpy(np.argsort(np.concatenate(bs))).to(feeder.dev)   # 길이순 -> smp 순
    return torch.cat(out)[back]
