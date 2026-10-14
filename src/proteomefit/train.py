"""훈련 루프 (분류 · 회귀 공통) · 손실 · 샘플 가중치"""
import numpy as np, pandas as pd, torch
from tqdm import tqdm
from .batching import make_batches


def taxon_weight(gen, org):
    """훈련 샘플의 genus · 균주 [m] -> 분류군 역크기 가중치 [m] = 1 / (genus 안 균주 수 x 균주 샘플 수)"""
    d = pd.DataFrame({'g': gen, 'o': org})
    return 1.0 / (d.groupby('g').o.transform('nunique') * d.groupby('o').o.transform('size')).to_numpy()


def ccc_loss(P, obs, w):
    """가지별 예측 P [B, K] · 라벨 obs [B] · 샘플 가중치 w [B] 또는 가지별 가중치 [B, K] (fp32) -> 가지별 가중 1 - CCC 의 평균 (스칼라)"""
    w = w if w.dim() == 2 else w[:, None]                          # [B, K] 또는 [B, 1]
    w = w / (w.sum(0) + 1e-8)                                      # 가지별 합 1 가중치
    o = obs[:, None]                                               # 라벨 [B, 1]
    mp, mo = (w * P).sum(0), (w * o).sum(0)                        # 가지별 가중 평균 [K] · 라벨 가중 평균 [K] (또는 [1])
    dp, do = P - mp, o - mo                                        # 편차 [B, K] · [B, K] (또는 [B, 1])
    vp, vo, cov = (w * dp * dp).sum(0), (w * do * do).sum(0), (w * dp * do).sum(0)   # 가중 분산 · 라벨 분산 · 공분산 [K]
    return (1 - 2 * cov / (vp + vo + (mp - mo) ** 2 + 1e-8)).mean()


def bootstrap_counts(tr, n, n_branch, seed, device):
    """훈련 샘플 번호 tr [m] -> 가지별 부트스트랩 등장 횟수 [n, n_branch] uint8 GPU (훈련 밖 샘플 = 0)
       가지마다 tr 에서 m 개를 복원 추출한다"""
    g = torch.Generator(device).manual_seed(seed)
    t = torch.from_numpy(tr).to(device)                            # 훈련 샘플 번호 [m]
    C = torch.zeros(n, n_branch, dtype=torch.uint8, device=device) # [n, n_branch]
    for k in range(n_branch):
        draw = torch.randint(len(t), (len(t),), generator=g, device=device)   # 복원 추출한 tr 안 위치 [m]
        C[t, k] = torch.bincount(draw, minlength=len(t)).to(torch.uint8)      # 등장 횟수 [m] (평균 1)
    return C


def train_model(task, inputs, feeder, lens, idx, seed, tag, cfg, sw=None, n_epoch=None):
    """입력 조합 inputs 의 모델 하나 훈련 -> (model, net, 히스토리 DataFrame, best epoch)
       lens: 배치 구성용 샘플 길이 [n] / idx: {'train', 'val'} 샘플 번호 / sw: 훈련 샘플 가중치 [n_train]
       n_epoch=None: 검증 지표(task.score)로 조기 종료하고 best epoch 가중치를 복원한다
       n_epoch 지정: 그만큼 훈련하고 마지막 가중치를 쓴다 (최종 모델)"""
    dev = cfg.device
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    model = task.build_model(inputs, seed).to(dev)
    net = torch.compile(model, dynamic=True) if cfg.compile else model    # model 과 파라미터 공유

    tr, ve = idx['train'], task.val_samples(idx['val'])            # 훈련 샘플 · 조기 종료 검증 샘플
    tr_score = task.scorer(tr)                                     # 훈련 지표 함수 (tr 순서 예측 -> 스칼라)
    early = n_epoch is None
    if early:
        vl_score = task.scorer(ve)
    sw_full = torch.ones(task.n, device=dev)                       # 전체 샘플 길이 가중치 [n]
    if sw is not None:
        sw_full[torch.from_numpy(tr).to(dev)] = torch.from_numpy(sw).to(dev)
    task.begin(tr, sw_full, seed)                                  # 손실 준비 (양성 가중 · 부트스트랩 횟수)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg.lr, weight_decay=cfg.wd, fused=True)

    sc = task.score
    hist, best_r, best_ep, best_state, pat = [], -2.0, 0, None, 0
    bar = tqdm(range(1, (n_epoch or cfg.max_epoch) + 1), desc=tag, mininterval=0.5)
    for ep in bar:
        model.train()
        bs = make_batches(tr, lens, cfg.batch, cfg.tok_budget, rng)
        run_loss, pts = torch.zeros((), device=dev), []
        for *xb, ds, b in feeder.stream(bs):                       # b = GPU 샘플 번호 [B]
            with torch.autocast('cuda', dtype=cfg.amp_dtype):
                out = net(*xb, ds)                                 # 로짓 [B] 또는 가지별 예측 [B, K]
            loss = task.batch_loss(out.float(), b)
            opt.zero_grad(set_to_none=True); loss.backward(); opt.step()
            run_loss += loss.detach() * len(b)
            pts.append(task.point(out.detach().float()))           # 샘플별 점 예측 [B]
        back = torch.from_numpy(np.argsort(np.concatenate(bs))).to(dev)   # 처리 순서 -> tr 순서
        tr_loss = (run_loss / len(tr)).item()
        tr_s    = tr_score(torch.cat(pts)[back])

        if early:
            raw = task.predict(net, feeder, lens, ve)              # 검증 예측
            vl_loss, vl_s = task.val_loss(raw, ve), vl_score(task.val_point(raw))
            if np.isfinite(vl_s) and vl_s > best_r + cfg.min_delta:
                best_r, best_ep, pat = vl_s, ep, 0
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                pat += 1
        else:
            vl_loss, vl_s = np.nan, np.nan
        hist.append({'epoch': ep, 'tr_loss': tr_loss, 'vl_loss': vl_loss, f'tr_{sc}': tr_s, f'vl_{sc}': vl_s})
        bar.set_postfix_str(f'loss={tr_loss:.4f}/{vl_loss:.4f} {sc}={tr_s:.3f}/{vl_s:.3f} '
                            f'best={best_r:.3f}@{best_ep} pat={pat}', refresh=False)
        if early and pat >= cfg.patience:
            break
    bar.close()
    if best_state is not None:
        model.load_state_dict(best_state)                          # best epoch 가중치 복원
    hist = pd.DataFrame(hist)
    hist['is_best'] = hist.epoch == best_ep
    return model, net, hist, best_ep

