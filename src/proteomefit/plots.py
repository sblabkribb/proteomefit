"""그림 - 학습 곡선 · 혼동행렬 · 실측-예측 산점도 · 신뢰도 보정 (파일로 저장하고 닫는다)"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix


def plot_history(hist, best_ep, path, score, loss_label, score_label):
    """epoch별 손실(왼쪽) · 조기 종료 지표(오른쪽, 열 tr_<score> · vl_<score>) 곡선, best epoch 점선"""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, (a, b, name) in zip(axes, [('tr_loss', 'vl_loss', loss_label),
                                       (f'tr_{score}', f'vl_{score}', score_label)]):
        ax.plot(hist.epoch, hist[a], '-o', lw=1.2, ms=2.5, label='train')
        ax.plot(hist.epoch, hist[b], '-o', lw=1.2, ms=2.5, label='val')
        if best_ep:
            ax.axvline(best_ep, ls='--', c='0.4', lw=1, label=f'best epoch {best_ep}')
        ax.set_xlabel('epoch'); ax.set_ylabel(name); ax.grid(alpha=.3)
    fig.tight_layout(rect=[0, .09, 1, 1])
    fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center', ncol=3, frameon=False)
    fig.savefig(path, dpi=140); plt.close(fig)


def confusion_plot(y_true, preds, path, title, classes):
    """분류 - 라벨 [m] · {패널 제목: 0/1 예측 [m]} -> 혼동행렬 [1, 패널 수] (행 = 예측 · 열 = 실제, 실제 클래스 기준 정규화)
       대제목 아래 줄 = 전체 n"""
    fig, axes = plt.subplots(1, len(preds), figsize=(5.2 * len(preds), 4.8), squeeze=False)
    for ax, (nm, pred) in zip(axes[0], preds.items()):
        C = confusion_matrix(pred, y_true, labels=[0, 1])                         # [2, 2] 예측(행) x 실제(열)
        R = C / np.maximum(C.sum(0, keepdims=True), 1)                          # 실제 클래스(열)별 비율
        ax.imshow(R, cmap='Blues', vmin=0, vmax=1)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f'{R[i, j]:.2f}\n{int(C[i, j]):,}', ha='center', va='center',
                        fontsize=9, color='white' if R[i, j] > .5 else 'black')
        ax.set_xticks(range(2), classes, fontsize=8); ax.set_yticks(range(2), classes, fontsize=8)
        ax.set_xlabel('observed'); ax.set_ylabel('predicted'); ax.set_title(nm, fontsize=9)
    fig.suptitle(f'{title}\nn = {len(y_true):,}', fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, .91]); fig.savefig(path, dpi=140); plt.close(fig)


def scatter_plot(sets, path, title):
    """회귀 - (이름, 실측, 예측 평균, 막대 반폭) 목록 -> 실측(가로) - 예측 평균(세로) 산점도 패널
       점 = 회색, 막대 = 예측 평균 ± 반폭, 검은 점선 = y = x, 빨간 점선 = 최소제곱 적합선. 축 범위는 첫 패널 기준 공통"""
    _, o0, p0, _ = sets[0]
    lim = max(np.abs(o0).max(), np.abs(p0).max()) * 1.08          # 공통 축 반폭 (첫 패널 기준)
    fig, axes = plt.subplots(1, len(sets), figsize=(6.2 * len(sets), 6.6))
    for ax, (nm, obs, pred, ci) in zip(axes, sets):
        a, b = np.polyfit(obs, pred, 1)                             # 적합선 기울기 · 절편
        ax.errorbar(obs, pred, yerr=ci, fmt='none', ecolor='#a0a0a0', elinewidth=.3, zorder=1)
        ax.scatter(obs, pred, s=6, c='#606060', linewidths=0, zorder=2)
        xs = np.array([-lim, lim])
        ax.plot(xs, xs, 'k--', lw=1.2)                              # y = x
        ax.plot(xs, a * xs + b, 'r--', lw=1.2)                      # 적합선
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect('equal')
        ax.set_xlabel('observed fitness (per-organism z)'); ax.set_ylabel('predicted mean (per-organism z)')
        ax.set_title(f'{nm}\nfit slope={a:.3f}  intercept={b:+.3f}', fontsize=9)
        ax.grid(alpha=.25)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, .95]); fig.savefig(path, dpi=140, bbox_inches='tight'); plt.close(fig)


def sign_confusion_plot(sets, path, title, classes):
    """회귀 - (이름, 라벨, 예측 평균) 목록 -> 0 기준 이진 분류 혼동행렬 (행 = 예측, 열 = 실제, 실제 클래스 기준 정규화)"""
    fig, axes = plt.subplots(1, len(sets), figsize=(4.6 * len(sets), 4.4))
    for ax, (nm, yt, mu) in zip(axes, sets):
        C = confusion_matrix(mu > 0, yt > 0, labels=[False, True])  # [2, 2] 예측(행) x 실제(열)
        R = C / C.sum(0, keepdims=True)                             # 실제 클래스(열)별 비율
        ax.imshow(R, cmap='Blues', vmin=0, vmax=1)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f'{R[i, j]:.2f}\n{C[i, j]:,}', ha='center', va='center', fontsize=9,
                        color='white' if R[i, j] > .5 else 'black')
        ax.set_xticks(range(2), classes); ax.set_yticks(range(2), classes)
        ax.set_xlabel('observed (y > 0)'); ax.set_ylabel('predicted (mean > 0)'); ax.set_title(nm, fontsize=9)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, .92]); fig.savefig(path, dpi=140); plt.close(fig)


def reliability_plot(tables, path, title='', notes=None, bins=np.linspace(0, 1, 11)):
    """{신뢰도 이름: metrics.reliability_table 결과} · 저장 경로 · 제목 · {신뢰도 이름: {단계: 지표 줄}} · 확률 구간 경계 -> 보정 그림
       신뢰도 1개 = [1, 2] (보정 전 | 후), 여러 개 = [2, 신뢰도 수] (위 행 = 보정 전, 아래 행 = 보정 후)
       구간마다 원 하나 - 가로 = 구간 예측 확률 평균, 세로 = 실제 비율 (관례), 회색 격자 = 구간 경계 (bins),
       원 넓이 = 구간 예측 수 / 전체 예측 수 (최소 4 pt^2). 대각선 아래 = 과대 예측, 위 = 과소 예측
       패널 제목 아래 줄 = notes 의 지표, 대제목 아래 줄 = 전체 n, 그림 아래 = 원 크기 범례 (n) 한 줄"""
    n = len(tables)
    N = next(iter(tables.values())).query("stage == 'after'").n.sum()        # 전체 예측 수 (패널 공통)
    size = lambda v: np.maximum(800 * v / N, 4)                              # 예측 수 -> 원 넓이 (pt^2)
    fig, axes = plt.subplots(*((1, 2) if n == 1 else (2, n)), figsize=((10, 5.9) if n == 1 else (4.8 * n, 10.3)), squeeze=False)
    for c, (name, T) in enumerate(tables.items()):
        for r, stage in enumerate(('before', 'after')):
            ax, t = (axes[0, r] if n == 1 else axes[r, c]), T[T.stage == stage]
            ax.scatter(t.predicted, t.observed, s=size(t.n), color=f'C{r}', edgecolor='white', zorder=3)   # C0 = 보정 전, C1 = 보정 후
            ax.plot([0, 1], [0, 1], 'k--', lw=1)
            ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')
            for e in bins:                                                             # 회색 격자 = 구간 경계 (가로 · 세로 모두, 눈금과 별개로 긋는다)
                ax.axvline(e, color='0.75', lw=0.6, zorder=1)
                ax.axhline(e, color='0.75', lw=0.6, zorder=1)                          # 눈금 숫자 = matplotlib 기본 (0.2 간격)
            ax.set_xlabel(f'predicted {name} (bin mean)'); ax.set_ylabel('observed share')
            note = (notes or {}).get(name, {}).get(stage, '')
            ax.set_title(f'{name} - {stage} calibration' + (f'\n{note}' if note else ''), fontsize=9)
    vals = [int(float(f'{N * q:.1g}')) for q in (0.01, 0.1, 0.5)]           # 범례 n = 전체의 약 1 · 10 · 50 % (유효숫자 1자리)
    hs = [plt.scatter([], [], s=size(v), color='0.55', edgecolor='white') for v in vals]
    fig.legend(hs, [f'{v:,}' for v in vals], title='n (predictions in the bin)', loc='lower center', ncol=len(vals),
               bbox_to_anchor=(0.5, 0.03), fontsize=8, title_fontsize=8, frameon=False, columnspacing=3, handletextpad=1)
    fig.text(0.5, 0.01, 'held-out genera; one circle per probability bin (grey grid lines = bin edges); '
             'circle area is proportional to the share of predictions in the bin',
             ha='center', fontsize=8, color='#555555')
    fig.suptitle((f'{title}\n' if title else '') + f'n = {N:,}', fontsize=10)
    fig.tight_layout(rect=(0, 0.15 if n == 1 else 0.08, 1, 0.93 if n == 1 else 0.96), h_pad=3); fig.savefig(path, dpi=140); plt.close(fig)
