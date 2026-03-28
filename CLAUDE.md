# 기본 설정

## 언어
- 항상 한국어로 답해줘
- 코드 주석은 한국어로 달아줘
- 에러 설명도 한국어로 해줘

## 코드 스타일
- 코드는 간결하고 읽기 쉽게 작성해줘
- 중요한 부분엔 주석을 달아줘

## 응답 스타일
- 설명은 단계별로 쉽게 해줘
- 모르는 용어는 풀어서 설명해줘

## 전공 및 코딩 스타일
- 나는 계량경제학(econometrics) 전공 대학원생이야
- 주 언어는 R이야
- 코드는 통계학자/계량경제학자가 작성한 것처럼 작성해줘
- 변수명은 논문 표기법을 따라줘 (예: y_hat, beta_ols, se_robust)
- 회귀분석 결과는 stargazer 또는 modelsummary로 정리해줘
- R 코드를 나에게 주기 전에 항상 Rscript로 직접 실행해서 에러 없이 돌아가는지 먼저 확인할 것

## 현재 연구 컨텍스트: Online Updating for Linear Panel Regressions

- 나는 현재 **Online Updating for Linear Panel Regressions**라는 계량경제학 방법론 연구를 진행 중이야.
- 연구 주제는 **선형 패널회귀(linear panel regression)**에서 새로운 데이터가 순차적으로 도착할 때, 전체 데이터를 다시 저장하거나 재추정하지 않고도 **계수 추정치와 표준오차를 정확하게 갱신하는 방법**을 개발하는 것이야.
- 이 연구의 핵심 문제의식은 다음과 같아:
  - 데이터가 매우 커서 전체 데이터를 메모리에 올릴 수 없거나
  - 보안/기밀성 때문에 과거 원자료 전체에 다시 접근할 수 없을 때
  - 이전에 저장한 **저차원 summary objects**만으로 추정량을 업데이트하고 싶다는 점이야.

### 연구 모형과 범위
- 주요 모형은 다음과 같아:
  - **Pooled OLS (POLS)**
  - **One-Way Fixed Effects (OWFE)**
  - **Two-Way Fixed Effects (TWFE)**
- 특히 패널데이터의 확장은 두 가지 경우를 구분해 다뤄:
  1. **새로운 individual이 들어오는 경우**
  2. **기존 individual에 새로운 time observation이 추가되는 경우**
- TWFE에서는 두 번째 경우를 다시 구분해야 해:
  - **기존 calendar time에 속한 새 관측치**인지
  - **완전히 새로운 calendar time (새 time dummy 필요)** 인지
- 이 차이를 매우 중요하게 다뤄야 해.

### 지금까지 정리된 핵심 결과
- OWFE와 TWFE에 대해 다음을 온라인으로 갱신하는 closed-form update를 유도했어:
  - **회귀계수**
  - **classical variance**
  - **cluster-robust variance**
- 알고리즘 구조는 다음처럼 이해하면 돼:
  - **Algorithm 1**: 새 individual 유입
  - **Algorithm 2**: 기존 individual에 기존 support 내 새 time observation 유입
  - **Algorithm 3**: 기존 individual에 새로운 calendar time 유입 (TWFE에서 차원 증가)
- 답변할 때는 가능하면 항상
  - 어떤 quantity를 저장해야 하는지
  - 무엇이 바뀌고 무엇은 그대로인지
  - 어떤 행렬/벡터를 업데이트해야 하는지
  - offline refit 없이 왜 가능한지
  를 분명히 설명해줘.

### 연구상 중요한 해석 포인트
- 이 연구는 **stochastic approximation**이나 **근사 알고리즘**이 아니라,
  **offline full-sample refit과 algebraically exact하게 일치하는 online updating**이 핵심이야.
- 특히 TWFE에서는 unbalanced panel 때문에 단순 double demeaning이 일반적으로 통하지 않아서,
  **individual FE는 within transformation으로 제거하고 time FE는 dummy regressor로 포함하는 방식**을 사용해.
- 새로운 calendar time이 들어오는 경우에는 time dummy가 하나 추가되어 **parameter dimension이 증가**한다는 점을 항상 주의해야 해.
- 이 경우 기존 객체를 단순히 확장할 수 있는 이유는,
  **새로 추가되는 time dummy가 pre-update sample에서는 항상 0**이기 때문이야.
- 반대로 임의의 새로운 regressor가 추가되는 상황과는 다르므로, 둘을 혼동하면 안 돼.

### 구현 및 검증 관련 맥락
- 주 구현 언어는 **R**이야.
- 검증은 주로 **plm** 패키지를 offline benchmark로 사용해 수행했어.
- 비교 대상은
  - coefficient
  - classical variance
  - cluster-robust variance
  이고,
  결과는 가능한 한 **machine precision 수준까지 일치하는지**를 확인하는 방식으로 본다.
- 따라서 코드를 작성하거나 수정할 때는
  - 수식과 동일한 notation 사용
  - 온라인 상태(state)에 저장해야 하는 객체 명시
  - batch refit과 비교 가능한 검증 코드
  를 중요하게 생각해줘.

### 답변 스타일에 대한 추가 지침
- 이 연구와 관련된 질문에는 일반적인 회귀 코드 답변보다
  **논문 notation 중심**으로 답해줘.
- 가능하면 다음 순서로 설명해줘:
  1. 문제 설정
  2. 업데이트 전 저장 상태(state)
  3. 새 데이터가 들어왔을 때 변하는 항
  4. coefficient update
  5. variance update
  6. 구현상 필요한 저장 객체
- 내가 논문 문장, 알고리즘 설명, 증명, R 코드, verification 코드를 요청하면
  항상 **계량경제학 논문 스타일**로 엄밀하게 답해줘.
- 수식 설명에서는 되도록 논문 표기법을 유지하고,
  막연한 설명보다 **정확히 어떤 행렬이 바뀌는지**를 중심으로 설명해줘.

### 특히 기억할 세부 선호
- OWFE 새 time observation update에서는
  **오직 해당 individual의 within transformation만 변한다**는 점을 강조해줘.
- classical variance update 설명에서는
  **pre-update inverse, sigma^2, coefficient, individual mean**이 왜 필요한지 분명히 설명해줘.
- cluster-robust variance update 설명에서는
  추가로 저장해야 하는 aggregate matrix가 무엇인지 구체적으로 적어줘.
- TWFE 새 calendar time update에서는
  **zero-extension**, **dimension expansion**, **새 time dummy의 pre-update 값이 모두 0**이라는 점을 반드시 언급해줘.
- 알고리즘이나 문장 수정 시에는
  “무엇을 업데이트해야 하는지”와 “왜 전체 데이터를 다시 볼 필요가 없는지”가 드러나게 써줘.
- 자세한 내용은 online_algorithms.pdf 논문 파일을 읽을 필요가 있어.