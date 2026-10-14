"""모델 - PMA 풀링 분기 + 연결 벡터 위의 MLP 머리(분류) 또는 미니 MLP 뱅크(회귀)

입력 조합 inputs 는 ('seq', 'esm', 'prot', 'bio', 'ctx', 'ctxp', 'ctxb') 같은 이름 목록이다.
seq · esm 은 가변 길이 입력(PMA 로 축약), 나머지는 데이터셋 블록이며 이 순서대로 이어 붙인다.
"""
import numpy as np, torch, torch.nn as nn

VAR_INPUTS = ('seq', 'esm')                                           # 가변 길이 입력


def ds_blocks(inputs):
    """입력 목록 -> 데이터셋 블록 목록 (esm · seq 제외, 순서 유지)"""
    return [b for b in inputs if b not in VAR_INPUTS]


class PMA(nn.Module):
    """Pooling by Multihead Attention (시드 1개) - 토큰 [B,T,D] + 패딩 마스크 [B,T] -> [B,D]
       파라미터 = nn.MultiheadAttention, 계산은 고정 시드 쿼리로 직접 한다"""
    def __init__(self, dim, n_head, droprate):
        super().__init__()
        self.q = nn.Parameter(torch.randn(1, 1, dim) * 0.02)          # 시드 벡터 [1,1,D]
        self.attn = nn.MultiheadAttention(dim, n_head, dropout=droprate, batch_first=True)
        self.norm = nn.LayerNorm(dim)

    def forward(self, tok, pad):
        """nn.MultiheadAttention(q, tok, tok) 와 같은 값 - 키 투영은 쿼리 쪽 [D, H] 가중치로 접고, 값 투영은 토큰 가중 평균 뒤에 한다"""
        B, D, H = len(tok), tok.shape[-1], self.attn.num_heads
        Wq, Wk, Wv = self.attn.in_proj_weight.view(3, H, D // H, D)          # 헤드별 투영 [H, dh, D] x 3
        bq, _, bv = self.attn.in_proj_bias.view(3, H, D // H)                # 헤드별 편향 [H, dh] x 3 (키 편향은 상쇄)
        qh = (torch.einsum('hkd,d->hk', Wq, self.q.view(D)) + bq) / (D // H) ** 0.5   # 스케일 포함 헤드별 쿼리 [H, dh]
        U  = torch.einsum('hk,hkd->dh', qh, Wk)                              # 키 투영을 쿼리 쪽으로 접은 가중치 [D, H]
        s  = (tok @ U.to(tok.dtype)).float().masked_fill(pad[..., None], float('-inf'))   # 점수 [B, T, H] (패딩 제외)
        a  = nn.functional.dropout(s.softmax(1), self.attn.dropout, self.training)   # 어텐션 가중치 [B, T, H]
        P  = torch.einsum('bth,btd->bhd', a.to(tok.dtype), tok)              # 헤드별 토큰 가중 평균 [B, H, D]
        o  = torch.einsum('bhd,hkd->bhk', P, Wv.to(P.dtype)) + a.sum(1)[..., None] * bv   # 값 투영 [B, H, dh]
        return self.norm(self.q.view(1, D) + self.attn.out_proj(o.reshape(B, D)))   # 잔차 + LayerNorm [B, D]


class SeqBranch(nn.Module):
    """정수 서열 [B,L] -> 임베딩 + 학습 위치 임베딩 -> PMA -> [B,dim]"""
    def __init__(self, dim, n_head, droprate, lmax, vocab=22):
        super().__init__()
        self.emb = nn.Embedding(vocab, dim, padding_idx=0)
        self.pos = nn.Parameter(torch.randn(1, lmax, dim) * 0.02)
        self.pool = PMA(dim, n_head, droprate)

    def forward(self, seq):
        return self.pool(self.emb(seq) + self.pos[:, :seq.shape[1]], seq == 0)


class Backbone(nn.Module):
    """입력 조합 inputs 의 연결 벡터 - seq · esm 분기(입력에 있을 때만) + 데이터셋 -> [B, d_in]
       block_dim: 입력별 폭 dict"""
    def __init__(self, inputs, block_dim, seq_dim=64, seq_head=4, n_head=8, droprate=0.3, max_len=4096):
        super().__init__()
        self.pool_esm = PMA(block_dim['esm'], n_head, droprate) if 'esm' in inputs else None
        self.seq = SeqBranch(seq_dim, seq_head, droprate, max_len) if 'seq' in inputs else None
        self.var = [b for b in inputs if b in VAR_INPUTS]             # 가변 길이 입력의 연결 순서
        self.widths = [block_dim[b] for b in self.var + ds_blocks(inputs)]   # 연결 순서의 블록 폭
        self.d_in = sum(self.widths)                                  # 연결 벡터 폭

    def concat(self, tok, pad, seq, ds):
        h = {}
        if self.pool_esm is not None:                                 # ESM 분기 먼저 (드롭아웃 난수 순서)
            h['esm'] = self.pool_esm(tok, pad).float()                # [B, 1152]
        if self.seq is not None:
            h['seq'] = self.seq(seq).float()                          # [B, 64]
        return torch.cat([h[b] for b in self.var] + [ds.float()], 1)  # 입력 목록 순서대로 연결 [B, d_in]


class EssNet(Backbone):
    """분류 모델 - 연결 벡터 -> Dropout -> [Linear-GELU-Dropout] x len(mlp_hidden) -> Linear(1) -> 로짓 [B] (sigmoid = P(essential))"""
    def __init__(self, inputs, block_dim, mlp_hidden=(640, 320), droprate=0.3, **kw):
        super().__init__(inputs, block_dim, droprate=droprate, **kw)
        d, layers = self.d_in, [nn.Dropout(droprate)]
        for h in mlp_hidden:
            layers += [nn.Linear(d, h), nn.GELU(), nn.Dropout(droprate)]
            d = h
        self.mlp = nn.Sequential(*layers, nn.Linear(d, 1))

    def forward(self, tok, pad, seq, ds):
        return self.mlp(self.concat(tok, pad, seq, ds)).squeeze(1)    # [B] 로짓


def block_bounds(widths):
    """블록 폭 목록 -> 블록 경계 [(start, stop), ...]"""
    e = np.cumsum([0] + list(widths))
    return [(int(a), int(b)) for a, b in zip(e[:-1], e[1:])]


def make_branch_idx(bounds, n_branch, m, rng, block_p):
    """블록 경계 [(start, stop), ...] -> 가지별 입력 열 번호 [n_branch, m] int64
       가지마다 블록을 확률 block_p 로 고르고(하나도 없으면 1개), 고른 블록에 열 m 개를 균등 배분해 무작위로 뽑는다"""
    idx = np.empty((n_branch, m), np.int64)
    for b in range(n_branch):
        pick = [k for k in range(len(bounds)) if rng.random() < block_p] or [int(rng.integers(len(bounds)))]
        share = np.full(len(pick), m // len(pick))
        share[:m - share.sum()] += 1                                  # 나머지 열을 앞 블록부터 1개씩
        cols = [rng.choice(np.arange(*bounds[k]), q, replace=bounds[k][1] - bounds[k][0] < q)
                for k, q in zip(pick, share)]                         # 블록이 할당보다 작으면 복원 추출
        idx[b] = np.concatenate(cols)
    return torch.from_numpy(idx)


class BranchBank(nn.Module):
    """미니 MLP n_branch 개를 배치 행렬곱으로 한 번에 계산 - 연결 벡터 [B, D] -> 가지별 예측 [B, n_branch]
       가지 = [Linear -> GELU -> Dropout] x len(hidden) -> Linear(1), 입력 열 = idx [n_branch, m]"""
    def __init__(self, idx, hidden, droprate):
        super().__init__()
        self.register_buffer('idx', idx)
        K, m = idx.shape
        self.drop = droprate
        dims = [m, *hidden, 1]                                        # 층별 노드 수
        self.W = nn.ParameterList([nn.Parameter(torch.randn(K, a, b) * ((2.0 if b > 1 else 1.0) / a) ** 0.5)
                                   for a, b in zip(dims[:-1], dims[1:])])   # 가지별 가중치 [K, a, b] (은닉층 He, 출력층 1/a)
        self.b = nn.ParameterList([nn.Parameter(torch.zeros(K, 1, b)) for b in dims[1:]])   # 가지별 편향 [K, 1, b]

    def forward(self, x):
        h = x[:, self.idx].permute(1, 0, 2)                           # 가지별 입력 [K, B, m]
        for W, b in zip(self.W[:-1], self.b[:-1]):                    # 은닉층 [K, B, hidden_i]
            h = nn.functional.dropout(nn.functional.gelu(torch.baddbmm(b, h, W)), self.drop, self.training)
        return torch.baddbmm(self.b[-1], h, self.W[-1]).squeeze(-1).T  # [B, K]


class FitNet(Backbone):
    """회귀 모델 - 연결 벡터 -> Dropout -> 미니 MLP 뱅크 -> 가지별 예측 [B, n_branch] (균주별 z-정규화 적합도)
       seed: 가지별 입력 열 선택 난수 시드"""
    def __init__(self, inputs, block_dim, seed, n_branch=1000, br_feat=128, br_hidden=(32,), block_p=0.5,
                 droprate=0.3, **kw):
        super().__init__(inputs, block_dim, droprate=droprate, **kw)
        self.drop = nn.Dropout(droprate)
        idx = make_branch_idx(block_bounds(self.widths), n_branch, br_feat, np.random.default_rng(seed), block_p)
        self.bank = BranchBank(idx, br_hidden, droprate)

    def forward(self, tok, pad, seq, ds):
        return self.bank(self.drop(self.concat(tok, pad, seq, ds)))   # [B, n_branch]


def param_summary(model):
    """모델 -> '전체 파라미터 수 = 자식 모듈별 파라미터 수' 문자열 (단위 백만)"""
    par = lambda m: sum(p.numel() for p in m.parameters()) / 1e6
    parts = [f'{nm} {par(m):.2f}M' for nm, m in model.named_children() if par(m) > 0]
    return f'{par(model):.2f}M = ' + ' + '.join(parts)
