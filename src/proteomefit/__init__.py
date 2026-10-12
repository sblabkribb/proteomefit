"""ProteomeFit - 단백질 언어 모델 임베딩 · 프로테옴 문맥으로 세균 유전자 필수성 · 적합도 예측

config   : 하이퍼파라미터 (Config) · torch 전역 설정
data     : 균주별 base data TSV · 임베딩 h5 적재
features : Biopython 서열 특성 · 정수 서열 · 직교상동체 목록 · 문맥 특성 · 데이터셋
folds    : genus 단위 중첩 교차검증 분할
models   : PMA · SeqBranch · EssNet(분류) · FitNet(회귀, 미니 MLP 뱅크)
batching : 배치 구성 · GPU 상주 배치 공급 · 추론
metrics  : 분류 · 회귀 지표
plots    : 학습 곡선 · 혼동행렬 · 산점도
train    : 훈련 루프 · 손실 · 샘플 가중치
tasks    : 과제 정의 (EssTask · FitTask)
cv       : 중첩 교차검증 실행 · 최종 모델 · 결과 요약
"""
