# BIOTRON Project — 생명 기원 시뮬레이션

## 프로젝트 소유자
권기현 (권앤컴퍼니 대표, 도호쿠대 화학과 졸업, AI 기업강사)

## 프로젝트 목표
**RL Agent가 PCR과 유사한 자기복제 프로토콜을 독자적으로 발견하는 것**

강화학습 에이전트가 원시 수프 시뮬레이션 환경에서 온도/에너지를 조절하여,
template-directed self-replication이 출현하는 조건(프로토콜)을 스스로 학습하게 한다.

## 배경 — 지금까지의 사고 과정

### 핵심 인사이트 (대화에서 도출)
1. **퍼셉트론 추상화**: 분자를 "같은 구조, 다른 값"으로 추상화. AGCU를 하드코딩하지 않음.
2. **D/A binding physics**: 결합 사이트는 Donor(1) 또는 Acceptor(0). D는 A에만 결합. 이것만으로 4종류 monomer 자동 출현 (α[00]↔ᾱ[11], β[01]↔β̄[10]).
3. **요철 모델**: Template가 monomer pool에서 complement를 "찍어냄". Template 자체가 촉매.
4. **Step-by-step assembly**: 한 번에 복제가 아님. Monomer가 하나씩 template에 결합. 옆에 이미 있으면 더 잘 붙음 (cooperative/stacking).
5. **에너지 사이클 = 자연의 PCR**: 고온에서 strand 분리(denaturation), 저온에서 monomer 결합(annealing). 이 사이클이 복제를 구동.
6. **볼츠만 생존**: P(break) = exp(-bondE/kT). 별도 decay 규칙 없이 열역학으로 통일.
7. **Birth sensitive, survival resilient**: 생명은 특정 조건에서만 탄생하지만, 한번 태어나면 다양한 조건에서 버팀. 환경 파라미터가 fluctuate해야 함.
8. **Template 효과는 코딩하지 않음**: H-bond + stacking + covalent polymerization 규칙만 넣으면 template-directed synthesis가 emergent하게 나와야 함.

### 5개의 물리 법칙 (이것만 시뮬레이션에 존재)
1. H-bond: D↔A (약, 가역)
2. Covalent bond: 인접 monomer 결합 (강, 비가역)
3. Stacking: 옆에 이미 붙어있으면 새 monomer가 더 잘 붙음
4. Energy: 환경 에너지 수준 (agent가 조절)
5. Boltzmann survival: P(break) = exp(-bondE/kT)

### 실험적 근거
- 밀러-유레이 실험: 아미노산, 뉴클레오베이스 합성 확인
- Joan Oró: HCN에서 adenine 합성
- Ferris: 점토 광물에서 50-mer RNA oligomer 합성
- Dimer/trimer 형성: prebiotic 조건에서 일상적으로 관찰
- Template-directed primer extension: trimer helper가 속도 100배 가속 (실험적 확인)
- Freeze-thaw cycle: cytosine, uracil 합성에 사용

## 현재 상태
- v1~v5d까지 React UI 시뮬레이터 개발 완료 (실험/시각화용)
- 아직 emergent self-replication 미달성
- 진단: 파라미터 탐색 공간이 너무 넓어서 사람이 수동으로 sweet spot을 찾을 수 없음

## 다음 단계: RL Agent

### Phase 1: Python 포팅
- biotron-v5d.jsx의 시뮬레이션 로직을 Python으로 포팅
- Headless (UI 없음), 빠른 실행
- 핵심 클래스: World, Strand, Monomer
- 5개 물리 법칙만 구현

### Phase 2: Gym 환경 래핑
```python
class BiotronEnv(gym.Env):
    # State: [strand_count, avg_length, max_length, hbond_count, 
    #         monomer_pool, current_energy, ...]
    # Action: energy_level (continuous 0~1, 또는 discrete [cold/warm/hot])
    # Reward: 
    #   +1 per birth
    #   +100 for sustained replication (N consecutive births)
    #   -0.01 per step (시간 비용)
```

### Phase 3: RL 학습
- Algorithm: PPO (stable-baselines3)
- Episodes: 100K+
- Episode length: 2000 steps
- Agent가 매 스텝 에너지 수준을 결정
- 목표: birth를 최대화하는 에너지 프로토콜 학습

### Phase 4: 분석
- Agent의 learned policy 분석
- 에너지 프로토콜 시각화: "이 순서로 온도를 바꾸면 복제가 일어난다"
- 실제 PCR 프로토콜과 비교
- 논문 초안: "RL Agent Independently Discovers PCR-like Protocol for Abiotic Self-Replication"

## 기술 스택
- Python 3.10+
- NumPy (시뮬레이션)
- Gymnasium (RL 환경)
- stable-baselines3 (PPO)
- matplotlib / plotly (분석)
- (optional) Ray RLlib for 대규모 학습

## 파일 구조
```
biotron-project/
├── CLAUDE.md          # 이 파일
├── sim/
│   ├── world.py       # 시뮬레이션 엔진
│   ├── strand.py      # Strand, Monomer 클래스
│   └── physics.py     # 5개 물리 법칙
├── env/
│   └── biotron_env.py # Gym 환경
├── train/
│   └── train_ppo.py   # RL 학습 스크립트
├── analysis/
│   └── analyze.py     # Policy 분석 및 시각화
├── react/
│   └── biotron-v5d.jsx # UI 시뮬레이터 (참고용)
└── README.md
```

## 핵심 원칙
- Template replication을 하드코딩하지 않는다 — 5개 물리 법칙에서 출현해야 한다
- AGCU를 하드코딩하지 않는다 — D/A physics에서 4문자가 출현해야 한다
- Agent는 에너지만 조절한다 — 나머지는 물리가 한다
- 발견된 프로토콜이 현실의 PCR과 유사하면 대성공
