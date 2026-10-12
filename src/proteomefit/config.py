"""실험 설정 - 03 · 04 노트북이 공유하는 하이퍼파라미터와 torch 전역 설정"""
from dataclasses import dataclass, field
import torch
import torch._dynamo, torch._inductor.config


@dataclass
class Config:
    # 경로 (노트북이 있는 src/ 기준)
    base_dir: str = '../data/basedata'               # 균주별 base data TSV 폴더 (01 산출물)
    prot_h5: str  = '../data/emb_proteome_328m.h5'   # 유전자별 ProteomeLM 임베딩 [n_all, 1152]
    mean_h5: str  = '../data/emb_esm_600m_mean.h5'   # 유전자별 ESM 토큰 평균 [n_all, 1152]
    esm_h5: str   = '../data/emb_esm_600m.h5'        # ESM 토큰별 임베딩 (유전자당 [T, 1152])
    out_dir: str  = '../results/esscls'              # 결과 폴더
    full_dir: str = ''                               # ablation 비교 기준(전체 입력) 결과 폴더

    # 라벨 · 평가
    label_col: str = 'fit_mean'                      # 라벨 원천 열
    t_eval: float  = 2.0                             # |comb_t| 기준 - 문맥 특성 otf_*strong* · 회귀 조기 종료 검증 샘플
    thresh: float  = 0.5                             # [분류] 양성 판정 임계 확률
    w_clip: tuple  = (0.3, 3.0)                      # [회귀] 훈련 샘플 가중치 clip(|comb_t|) 범위
    pt_q: float    = 0.10                            # [회귀] 주 지표 _pt 의 예측 t 상위 비율

    # 모델
    max_len: int    = 4096                           # 정수 서열 잔기 상한
    seq_dim: int    = 64                             # 정수 서열 분기 임베딩 차원
    seq_head: int   = 4                              # 정수 서열 분기 PMA 헤드 수
    n_head: int     = 8                              # ESM 분기 PMA 헤드 수 (1152의 약수)
    droprate: float = 0.3                            # 연결 벡터 · 은닉층 · PMA 드롭아웃
    mlp_hidden: tuple = (640, 320)                   # [분류] MLP 은닉 유닛 수
    n_branch: int   = 1000                           # [회귀] 미니 MLP 수
    br_feat: int    = 128                            # [회귀] 가지당 입력 열 수
    br_hidden: tuple = (32,)                         # [회귀] 미니 MLP 은닉층별 노드 수
    block_p: float  = 0.5                            # [회귀] 가지가 각 입력 블록을 고를 확률
    bootstrap: bool = True                           # [회귀] 가지마다 훈련 샘플 부트스트랩

    # 문맥 특성
    topk: tuple    = (1, 3, 5, 10, 20)               # 참조 상위 개수
    sim_chunk: int = 512                             # 직교상동체 검색의 질의 묶음 크기

    # 학습
    lr: float = 1e-3                                 # 학습률 (AdamW)
    wd: float = 1e-2                                 # 가중치 감쇄 (AdamW)
    batch: int = 1024                                # 훈련 배치 최대 샘플 수
    eval_batch: int = 4096                           # 평가 배치 최대 샘플 수
    tok_budget: int = 524288                         # 훈련 배치의 (샘플 수 x 최대 토큰 수) 상한
    eval_tok_budget: int = 1048576                   # 평가 배치의 같은 상한
    inner_k: int = 3                                 # 내부 폴드 수 = 검증 genus 수
    max_epoch: int = 500                             # 최대 epoch
    patience: int = 20                               # 조기 종료 인내 epoch
    min_delta: float = 0.01                          # 조기 종료 지표 개선 인정 하한
    seed: int = 42                                   # 초기화 시드
    csv_round: int = 6                               # csv 실수 반올림 자릿수
    resume: bool = True                              # 산출물이 완비된 폴드 · 최종 모델은 건너뛴다

    # 실행 환경
    device: str = field(default_factory=lambda: 'cuda:0' if torch.cuda.is_available() else 'cpu')
    amp_dtype: torch.dtype = torch.bfloat16          # 혼합정밀도 dtype
    t_align: int = 8                                 # 배치 내 최대 길이를 이 배수로 올림
    compile: bool = True                             # torch.compile 사용
    esm_chunk: int = 2048                            # ESM h5 를 한 번에 읽는 단백질 수
    gpu_reserve: int = 16                            # ESM 적재 후 남겨 둘 최소 VRAM (GiB)


def setup_torch():
    """torch 전역 설정 - TF32 행렬곱 · 동적 모양 컴파일 오류 회피"""
    torch.set_float32_matmul_precision('high')                        # autocast 밖 fp32 행렬곱을 TF32 로
    torch._dynamo.config.cache_size_limit = 64                        # 배치 모양별 컴파일 그래프 수 상한
    torch._inductor.config.triton.coalesce_tiling_analysis = False    # 동적 모양에서 sympy 검사 오류를 내는 분석 끔
    torch._inductor.config.shape_padding = False                      # 패딩 보폭 고정으로 인한 역전파 assert_size_stride 오류 회피
