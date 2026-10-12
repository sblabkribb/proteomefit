"""그림 - 학습 곡선 · 혼동행렬 · 실측-예측 산점도 (파일로 저장하고 닫는다)"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix


def plot_history(hist, best_ep, path, score, loss_label, score_label):
    """epoch별 손실(왼쪽) · 조기 종료 지표(오른쪽, 열 tr_<score> · vl_<score>) 곡선, best epoch 점선"""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, (a, b, name) in zip(axes, [('tr_loss', 'vl_loss', loss_label),
                                       (f'tr_{score}', f'vl_{score}', score_label)]):
        ax.plot(hist.epoch, hist[a], '-o', lw=1.2, ms=2.5, label='train')   # 마커: epoch 가 1개여도 보인다
        ax.plot(hist.epoch, hist[b], '-o', lw=1.2, ms=2.5, label='val')
        if best_ep:
            ax.axvline(best_ep, ls='--', c='0.4', lw=1, label=f'best epoch {best_ep}')
        ax.set_xlabel('epoch'); ax.set_ylabel(name); ax.grid(alpha=.3)
    fig.tight_layout(rect=[0, .09, 1, 1])
    fig.legend(*axes[0].get_legend_handles_labels(), loc='lower center', ncol=3, frameon=False)
    fig.savefig(path, dpi=140); plt.close(fig)


def ess_plot_title(head, n, df):
    """분류 그림 제목 두 줄 - df(균주별 지표) 지표 열의 평균"""
    return (f'{head}  |  n={n:,}\n'
            f'MCC {df.mcc.mean():.3f}  F1 {df.f1.mean():.3f}  recall {df.recall.mean():.3f}  '
            f'precision {df.precision.mean():.3f}  ROC AUC {df.roc_auc.mean():.3f}  PR AUC {df.pr_auc.mean():.3f}')


def prob_confusion_plot(y_true, prob, path, title, classes, thresh):
    """분류 - 혼동행렬(왼쪽, 행 = 예측 · 열 = 실제, 실제 클래스 기준 정규화) · 실제 클래스별 예측 확률 분포(오른쪽)"""
    C = confusion_matrix((prob >= thresh).astype(int), y_true, labels=[0, 1])   # [2, 2] 예측(행) x 실제(열)
    R = C / np.maximum(C.sum(0, keepdims=True), 1)                              # 실제 클래스(열)별 비율 (대각 = 클래스별 재현율)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    ax = axes[0]
    ax.imshow(R, cmap='Blues', vmin=0, vmax=1)
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f'{R[i, j]:.2f}\n{int(C[i, j]):,}', ha='center', va='center',
                    fontsize=8, color='white' if R[i, j] > .5 else 'black')
    ax.set_xticks(range(2), classes, fontsize=8); ax.set_yticks(range(2), classes, fontsize=8)
    ax.set_xlabel('observed'); ax.set_ylabel('predicted')
    ax.set_title('confusion (normalised by observed class)', fontsize=9)
    ax = axes[1]
    for c, nm in enumerate(classes):
        m = y_true == c
        if m.sum():
            ax.hist(prob[m], bins=np.linspace(0, 1, 51), density=True, histtype='step', lw=1.3,
                    label=f'{nm} (n={int(m.sum()):,})')
    ax.axvline(thresh, ls='--', c='0.4', lw=1)
    ax.set_xlabel(f'predicted P({classes[1]})'); ax.set_ylabel('density')
    ax.set_xlim(0, 1); ax.grid(alpha=.25); ax.legend(fontsize=7, framealpha=.8)
    ax.set_title('probability by observed class', fontsize=9)
    fig.suptitle(title, fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, .93])
    fig.savefig(path, dpi=140); plt.close(fig)


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
