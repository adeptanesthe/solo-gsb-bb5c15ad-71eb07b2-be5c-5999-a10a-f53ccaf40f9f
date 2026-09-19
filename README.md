# MSPC — 质谱同位素分布预计算 API

基于 **Python + FastAPI + SQLite** 的质谱预计算服务。客户端提交含元素计数的分子式（可带同位素标记），服务按各元素天然同位素丰度将分布表示为"质量-强度多项式"，逐元素卷积得到完整同位素包络，合并质量差小于 0.0001 Da 的峰，按基峰归一化后返回前 N 个峰及分子单同位素质量。

## 快速启动（单条命令）

```bash
docker compose up --build
```

启动过程会自动完成依赖安装（镜像构建时）与 SQLite 数据库初始化（服务启动时建表），无需任何手工迁移。

## 访问地址

| 地址 | 说明 |
|---|---|
| http://localhost:8000 | 服务根地址（端口 8000） |
| http://localhost:8000/docs | Swagger 交互式文档 |
| http://localhost:8000/health | 健康检查 |

## API 用法

### `POST /api/v1/envelope`

请求体：

```json
{
  "formula": "C6H12O6",
  "top_n": 5
}
```

- `formula`（必填）：分子式，如 `C6H12O6`；同位素标记用方括号表示，如 `C6H6[13C]2`（2 个纯 ¹³C 原子）、`C2H5[15N]O[18O]`。
- `top_n`（可选，默认 20，最大 200）：返回强度最高的前 N 个峰。

响应：

```json
{
  "formula": "C6H12O6",
  "composition": {"C": 6, "H": 12, "O": 6},
  "labels": {},
  "monoisotopic_mass": 180.063388,
  "base_peak_intensity": 100.0,
  "peak_count": 5,
  "peaks": [
    {"mass": 180.063388, "intensity": 100.0},
    {"mass": 181.066743, "intensity": 6.489437}
  ],
  "cached": false
}
```

- `monoisotopic_mass`：单同位素质量（各元素取天然丰度最高的同位素，标记原子取标记同位素质量）。
- `peaks`：基峰强度归一化为 100，按质量升序返回前 N 个峰。
- `cached`：结果是否来自 SQLite 缓存（相同 `formula` + `top_n` 的重复请求直接命中缓存）。

### `POST /api/v1/ion-envelope`

在中性分子式上施加**加合物**并返回带电离子的 m/z 同位素包络：先按加合物片段增减原子得到离子组成，复用现有分子式解析与多项式卷积生成中性质量包络，再转换为 m/z。

请求体：

```json
{
  "formula": "C6H12O6",
  "adduct": "+H",
  "charge": 1,
  "top_n": 5
}
```

- `formula`（必填）：中性分子式，规则同 `/api/v1/envelope`（可含 `[13C]` 等纯同位素标记）。
- `adduct`（必填）：由连续的"+/-[系数]分子式片段"组成，每个片段必须以显式 `+` 或 `-` 开头；系数为正整数（缺省为 1），**乘整个片段**。片段本身的语法与分子式相同（元素计数、`[13C]` 标记等），整个字符串必须被完整消费。例如：
  - `+H`：质子化（[M+H]⁺）；
  - `+Na-H`：钠离子加合物减一个氢；
  - `+2H-H2O`：加两个氢并脱去一个水分子（系数 2 作用于整个 `H`，即加 2 个 H）；
  - `-[13C]`：减去一个纯 ¹³C 原子。
  - 应用片段后任一元素（或同位素标记）计数变为负数时返回 422（零允许），`position` 指向导致**首次**负计数的片段首字符（`+`/`-` 与系数之后的第一个字符）。
- `charge`（必填）：带符号整数电荷，取值 **-10～-1 或 1～10**（不允许 0）。
- `top_n`（可选，默认 20，最大 200）：返回强度最高的前 N 个峰。

计算流程：

1. 用 `|charge| × MSPC_MERGE_TOLERANCE` 作为中性质量峰的合并容差（卷积及最终合并均用该放大容差，按默认 0.0001 Da 即电荷 2 时为 0.0002 Da）。
2. 对调整后的离子组成复用同一套卷积/剪枝/归一化（基峰 = 100）。
3. 每个峰按下式从离子质量（Da）转为 m/z（正离子扣除电子质量，负离子加回电子质量）：

   ```
   m/z = (mass - charge × 0.000548579909065 Da) / |charge|
   ```

响应：

```json
{
  "formula": "C6H12O6",
  "adduct": "+H",
  "charge": 1,
  "composition": {"C": 6, "H": 13, "O": 6},
  "labels": {},
  "monoisotopic_mz": 181.070665,
  "base_peak_intensity": 100.0,
  "peak_count": 5,
  "peaks": [
    {"mz": 181.070665, "intensity": 100.0},
    {"mz": 182.07402, "intensity": 6.489437}
  ]
}
```

- `composition` / `labels`：应用加合物并乘系数之后的离子组成（规则同 `/envelope` 的 `labels` 序列化）。
- `monoisotopic_mz`：调整后组成的单同位素中性质量经上式转换得到的单同位素 m/z。
- `peaks`：基峰归一化为 100 的前 N 个 m/z 峰，按 m/z 升序；双键电荷下相邻同位素峰的 m/z 间距约为 1/|charge| Da。

示例：`C6H12O6` + `+2H`、`charge` 2 的单同位素 m/z 为 `(180.063388 + 2×1.007825 - 2×0.0005486) / 2 ≈ 91.038971`；`-H`、`charge` -1 的去质子负离子 m/z 约为 `179.056112`。

服务默认监听 **8000 端口**（见 [访问地址](#访问地址)，容器内外均为 `http://localhost:8000`），直接用 curl 请求示例：

```bash
curl -X POST http://localhost:8000/api/v1/ion-envelope \
  -H 'Content-Type: application/json' \
  -d '{"formula":"C6H12O6","adduct":"+Na-H","charge":1,"top_n":5}'
```

### `POST /api/v1/ion-identify`

在一个中性分子式上同时评估 **1～20 个候选加合物/电荷组合**：对每个候选复用 `/api/v1/ion-envelope` 的加合物解析与离子包络计算（含 `|charge| ×` 合并容差缩放与 m/z 换算、基峰归一化为 100），将预测 m/z 峰与观测 m/z 峰在容差内一对一配对，以全部配对（含未匹配峰）的强度均方误差对候选排序，返回全部候选的得分、配对明细与最佳候选。

请求体：

```json
{
  "formula": "C6H12O6",
  "peaks": [
    {"mz": 181.070665, "intensity": 100.0},
    {"mz": 182.07402, "intensity": 6.49}
  ],
  "mass_tolerance": 0.01,
  "candidates": [
    {"adduct": "+Na-H", "charge": 1},
    {"adduct": "+H", "charge": 1},
    {"adduct": "+H", "charge": 2},
    {"adduct": "-99H", "charge": -1}
  ]
}
```

- `formula`（必填）：中性分子式，规则同 `/api/v1/ion-envelope`（可含 `[13C]` 等纯同位素标记）。
- `peaks`（必填，**1～200 个**）：观测峰，字段为 `mz`/`intensity`；m/z 必须**严格递增**，强度非负且至少一个为正。观测峰在评分前按基峰归一化为 100。
- `mass_tolerance`（必填，>0）：观测峰与预测 m/z 峰配对的质量容差（m/z 单位）。
- `candidates`（必填，**1～20 个**）：候选加合物/电荷组合，每个包含：
  - `adduct`（必填）：加合物字符串，语法与 `/api/v1/ion-envelope` 完全一致（`+/-[系数]片段`，如 `+H`、`+2H-H2O`、`+Na-H`）。
  - `charge`（必填）：带符号电荷，取值 **-10～-1 或 1～10**（不允许 0）。

配对与评分规则（与富集拟合接口同一套实现）：所有容差内的（观测, 预测）候选对按质量差绝对值从小到大贪心占用，每个峰至多参与一对；未匹配的观测峰以预测强度 0 计，未匹配的预测峰以观测强度 0 计；MSE 为全部配对强度差平方的均值。候选预测包络不经 `top_n` 截断（卷积剪枝阈值仍然生效）。

响应：

```json
{
  "formula": "C6H12O6",
  "mass_tolerance": 0.01,
  "observed_peak_count": 2,
  "best_candidate": {"index": 1, "adduct": "+H", "charge": 1, "mse": 0.0003},
  "candidates": [
    {
      "index": 1,
      "valid": true,
      "adduct": "+H",
      "charge": 1,
      "mse": 0.0003,
      "monoisotopic_mz": 181.070665,
      "composition": {"C": 6, "H": 13, "O": 6},
      "labels": {},
      "pair_count": 103,
      "pairs": [
        {"observed_mz": 181.070665, "observed_intensity": 100.0, "predicted_mz": 181.070665, "predicted_intensity": 100.0}
      ]
    },
    {
      "index": 0,
      "valid": true,
      "adduct": "+Na-H",
      "charge": 1,
      "mse": 185.99,
      "monoisotopic_mz": 202.044784,
      "composition": {"C": 6, "H": 11, "O": 6, "Na": 1},
      "labels": {},
      "pair_count": 108,
      "pairs": []
    },
    {
      "index": 3,
      "valid": false,
      "adduct": "-99H",
      "charge": -1,
      "mse": null,
      "monoisotopic_mz": null,
      "composition": null,
      "labels": null,
      "pair_count": null,
      "pairs": null,
      "error": "invalid_adduct",
      "message": "Adduct fragment would drive element 'H' count from 12 to -87; counts must not become negative",
      "position": 3
    }
  ]
}
```

- `candidates`：**全部候选**的结果，顺序为有效候选按 `mse` 升序、MSE 相同时按请求索引（输入顺序）稳定排序；非法候选排在末尾，彼此仍按请求索引升序。
  - 有效候选：`valid: true`，含 `mse`、`monoisotopic_mz`、`composition`、`labels`（含义同 `/api/v1/ion-envelope`）以及 `pairs` 配对明细。`pairs` 按 m/z 升序（未匹配预测峰按预测 m/z 插入）；未匹配观测峰的 `predicted_*` 为 `null`，未匹配预测峰的 `observed_*` 为 `null`。`index` 为该候选在请求中的 **0 基索引**（排序后保持不变）。
  - 非法候选：`valid: false`，评分相关字段为 `null`，并带 `error`（`invalid_adduct`）、`message` 与 `position`——即加合物解析既有的位置错误（含缺少 `+`/`-` 前导符、空白、片段内分子式错误，以及应用片段后计数变负，`position` 语义与 `/api/v1/ion-envelope` 完全一致）。非法候选**不会**导致整个请求 422。
- `best_candidate`：MSE 最小的候选（平局取输入顺序最前者），给出其请求 `index`、`adduct`、`charge`、`mse`；**全部候选均非法时为 `null`**。
- 请求级非法输入（分子式非法、峰列表为空/非递增/全零强度、`mass_tolerance` 非正、候选数不在 1～20、候选电荷为 0 或超出 ±10、峰数超过 200 等）仍按既有规则返回 422。

```bash
curl -X POST http://localhost:8000/api/v1/ion-identify \
  -H 'Content-Type: application/json' \
  -d '{"formula":"C6H12O6","mass_tolerance":0.01,
       "peaks":[{"mz":181.0707,"intensity":100.0},{"mz":182.074,"intensity":6.5}],
       "candidates":[{"adduct":"+H","charge":1},{"adduct":"+Na-H","charge":1}]}'
```

### `POST /api/v1/ion-mixture-fit`

在一个中性分子式上枚举 **2～10 个候选加合物/电荷** 的两两混合：每个有效候选复用 `/api/v1/ion-envelope` 的离子包络计算并各自按基峰归一化（基峰 = 100）；对每一对不同的有效候选，在首候选占比 0%～100% 的比例网格上构造混合谱，将混合谱与观测 m/z 峰用与 `/api/v1/ion-identify` 相同的一对一配对（含双侧未匹配峰）评分，返回最佳组合、全部有效组合排名、比例、MSE 与最佳配对明细。

混合谱的构造（对首候选占比 p、次候选占比 100%−p）：

1. 两个包络**保持各自的基峰归一化**（基峰 = 100），强度分别乘以 p/100 与 (100−p)/100 后叠加；
2. 叠加后质量差小于 `mass_tolerance` 的峰用既有峰合并规则融合（质量取强度加权平均、强度求和）；
3. 融合后的混合谱**再次按基峰归一化**（基峰 = 100）。端点 p=0%/100% 时占比为 0 的一侧不贡献峰，混合谱退化为另一候选的单包络；
4. 混合谱与同样按基峰归一化的观测峰做容差内一对一贪心配对，未匹配的观测峰以预测强度 0 计、未匹配的预测峰以观测强度 0 计，MSE 为全部配对强度差平方的均值。候选包络不经 `top_n` 截断（卷积剪枝阈值仍然生效）。

请求体：

```json
{
  "formula": "C6H12O6",
  "peaks": [
    {"mz": 181.070665, "intensity": 100.0},
    {"mz": 182.07402, "intensity": 6.49},
    {"mz": 203.0526, "intensity": 28.0}
  ],
  "mass_tolerance": 0.01,
  "candidates": [
    {"adduct": "+H", "charge": 1},
    {"adduct": "+Na-H", "charge": 1},
    {"adduct": "+2H", "charge": 2}
  ],
  "step_percent": 5.0
}
```

- `formula`（必填）：中性分子式，规则同 `/api/v1/ion-envelope`（可含 `[13C]` 等纯同位素标记）。
- `peaks`（必填，**1～200 个**）：观测 m/z 峰，规则同 `/api/v1/ion-identify`（严格递增、强度非负且至少一个为正），评分前按基峰归一化为 100。
- `mass_tolerance`（必填，>0）：m/z 容差，**同时**用于混合谱叠加时的峰融合与混合谱/观测峰的配对。
- `candidates`（必填，**2～10 个**）：候选加合物/电荷组合，字段与校验（`adduct` 语法、`charge` 为 ±1～±10）同 `/api/v1/ion-identify`。枚举范围为每两个**不同的有效候选**构成的无序对，即 n 个有效候选对应 C(n,2) 个组合。
- `step_percent`（可选，默认 1.0，范围 0.5～10）：首候选混合比例的网格步长（百分点），候选点为 0、step、2·step、…、100（100% 端点始终包含）。

平局规则：同一候选对内不同比例的 MSE 相同（浮点噪声以内）时取**较低的首候选比例**（网格按升序访问）；不同候选对之间 MSE 相同时取**较小的候选索引对**（按 `(首候选索引, 次候选索引)` 字典序，索引对始终满足首索引 < 次索引）。

响应：

```json
{
  "formula": "C6H12O6",
  "mass_tolerance": 0.01,
  "step_percent": 5.0,
  "observed_peak_count": 3,
  "valid_candidate_count": 3,
  "combination_count": 3,
  "grid_points": 21,
  "search_points": 63,
  "best_mixture": {
    "first_index": 0,
    "second_index": 1,
    "first_adduct": "+H",
    "first_charge": 1,
    "second_adduct": "+Na-H",
    "second_charge": 1,
    "first_fraction_percent": 70.0,
    "mse": 12.34
  },
  "combinations": [
    {
      "rank": 1,
      "first_index": 0,
      "second_index": 1,
      "first_adduct": "+H",
      "first_charge": 1,
      "second_adduct": "+Na-H",
      "second_charge": 1,
      "first_fraction_percent": 70.0,
      "mse": 12.34,
      "pair_count": 6,
      "pairs": [
        {"observed_mz": 181.070665, "observed_intensity": 100.0, "predicted_mz": 181.070665, "predicted_intensity": 100.0},
        {"observed_mz": null, "observed_intensity": null, "predicted_mz": 91.038971, "predicted_intensity": 3.1}
      ]
    }
  ],
  "candidates": [
    {"index": 0, "valid": true, "adduct": "+H", "charge": 1, "error": null, "message": null, "position": null},
    {"index": 1, "valid": false, "adduct": "-99H", "charge": -1, "error": "invalid_adduct", "message": "Adduct fragment would drive element 'H' count from 12 to -87; counts must not become negative", "position": 3}
  ]
}
```

- `best_mixture`：全部有效候选对中 MSE 最小的（候选对, 首候选比例）组合，给出两个候选的请求索引 `first_index`/`second_index`（前者较小）、加合物与电荷、`first_fraction_percent`（首候选占比百分比）与 `mse`；**有效候选少于两个时为 `null`**。
- `combinations`：**全部有效候选对**的排名，按 MSE 升序、平局按较小索引对稳定排序；`rank` 从 1 开始。每项含与 `best_mixture` 相同的组合/比例/MSE 字段，以及该组合最佳比例下的 `pair_count` 与 `pairs` 配对明细（字段与排序规则同 `/api/v1/ion-identify`：按 m/z 升序插入未匹配预测峰，未匹配观测峰的 `predicted_*` 为 `null`，未匹配预测峰的 `observed_*` 为 `null`）。
- `candidates`：**全部请求候选**的合法性结果（保持请求顺序），字段为 `index`/`valid`/`adduct`/`charge` 及非法时的 `error`（`invalid_adduct`）、`message`、`position`——与 `/api/v1/ion-identify` 完全相同的位置错误语义。**单个候选非法时返回其错误并跳过所有涉及该候选的组合**，不影响其他候选对评分，也不会导致整个请求 422。
- `valid_candidate_count`：有效候选数；`combination_count` 为有效候选对数 C(n,2)；`grid_points` 为每对尝试的比例点数；`search_points = combination_count × grid_points` 为总评分点数。
- 请求级非法输入（分子式非法、峰列表为空/非递增/全零强度、`mass_tolerance` 非正、候选数不在 2～10、候选电荷为 0 或超出 ±10、`step_percent` 越界、峰数超过 200 等）按既有规则返回 422。

```bash
curl -X POST http://localhost:8000/api/v1/ion-mixture-fit \
  -H 'Content-Type: application/json' \
  -d '{"formula":"C6H12O6","mass_tolerance":0.01,"step_percent":5.0,
       "peaks":[{"mz":181.0707,"intensity":100.0},{"mz":203.0526,"intensity":30.0}],
       "candidates":[{"adduct":"+H","charge":1},{"adduct":"+Na-H","charge":1}]}'
```

### `POST /api/v1/enrichment-fit`

估计分子中单个同位素标记的丰度（富集度）：在 0%～100% 范围内按给定步长网格搜索候选丰度，对每个候选值复用卷积生成归一化包络，与观测峰配对比较，返回均方误差最小的丰度。

请求体：

```json
{
  "formula": "C6",
  "label": "13C",
  "peaks": [
    {"mass": 72.0, "intensity": 5.0},
    {"mass": 73.003355, "intensity": 30.0},
    {"mass": 74.00671, "intensity": 75.0},
    {"mass": 75.010065, "intensity": 100.0},
    {"mass": 76.013419, "intensity": 75.0},
    {"mass": 77.016774, "intensity": 30.0},
    {"mass": 78.020129, "intensity": 5.0}
  ],
  "step_percent": 0.5,
  "mass_tolerance": 0.05
}
```

- `formula`（必填）：分子式，规则同 `/api/v1/envelope`（可含 `[13C]` 等纯同位素标记，标记贡献固定质量位移）。
- `label`（必填）：待估计的单个同位素标签，格式为"质量数+元素符号"，如 `13C`、`15N`、`2H`。
- `peaks`（必填，1～200 个）：观测峰，质量必须**严格递增**，强度非负且至少一个为正。
- `step_percent`（可选，默认 1.0，范围 0.5～10）：丰度搜索步长（百分点）。候选点为 0、step、2·step、…、100（100% 端点始终包含）。
- `mass_tolerance`（必填，>0）：观测峰与预测峰配对的质量容差（Da）。

拟合模型：候选丰度为 p 时，被标记元素的同位素分布改为"标记同位素占 p，其余天然同位素按比例分享 1−p"（p 等于天然丰度时即还原为天然分布），再与其他元素的天然分布卷积得到预测包络；预测包络与观测峰均按各自基峰归一化为 100 后比较。

响应：

```json
{
  "formula": "C6",
  "label": "13C",
  "best_abundance_percent": 50.0,
  "mse": 0.0,
  "search_points": 201,
  "pair_count": 7,
  "pairs": [
    {"observed_mass": 72.0, "observed_intensity": 5.0, "predicted_mass": 72.0, "predicted_intensity": 5.0},
    {"observed_mass": 73.003355, "observed_intensity": 30.0, "predicted_mass": 73.003355, "predicted_intensity": 30.0},
    {"observed_mass": 74.00671, "observed_intensity": 75.0, "predicted_mass": 74.00671, "predicted_intensity": 75.0},
    {"observed_mass": 75.010065, "observed_intensity": 100.0, "predicted_mass": 75.010065, "predicted_intensity": 100.0}
  ]
}
```

- `best_abundance_percent`：最佳丰度（百分比）。MSE 相同时取较低丰度。
- `mse`：最佳丰度下全部配对的强度均方误差。
- `search_points`：搜索点数（0%～100% 按步长划分的候选数）。
- `pairs`：最佳丰度下的逐峰配对，按质量升序。每个观测峰与容差内最近的未占用预测峰配对；未匹配的观测峰 `predicted_*` 为 `null`（预测强度按 0 计），未匹配的预测峰 `observed_*` 为 `null`（观测强度按 0 计）。两侧强度均为基峰=100 的归一化值。

配对与误差细节：所有容差内的（观测, 预测）候选对按质量差从小到大贪心占用，每个峰至多参与一对；MSE 在全部配对（含未匹配峰）上计算。

### `POST /api/v1/dual-enrichment-fit`

联合估计同一分子中**两个不同元素**的目标同位素丰度：对两个标签各自在 0%～100% 范围内按各自步长网格搜索，构成二维候选网格；每个候选点组合两种富集分布、复用卷积生成归一化包络，再用与 `/api/v1/enrichment-fit` 完全相同的一对一峰配对与含未匹配峰的 MSE 评分。

请求体：

```json
{
  "formula": "C2H5NO2",
  "labels": [
    {"label": "13C", "step_percent": 1.0},
    {"label": "15N", "step_percent": 1.0}
  ],
  "peaks": [
    {"mass": 75.032028, "intensity": 42.86},
    {"mass": 76.035383, "intensity": 100.0}
  ],
  "mass_tolerance": 0.05
}
```

- `formula`（必填）：分子式，规则同 `/api/v1/enrichment-fit`。
- `labels`（必填，**恰好两个**）：两个待估计同位素标签，每个包含：
  - `label`（必填）：同位素标签，格式同单标签接口，如 `13C`、`15N`。
  - `step_percent`（可选，默认 1.0，范围 0.5～10）：该标签自己的丰度搜索步长。
- 两个标签必须属于**不同元素**（如 `13C` 与 `12C` 同属 C，属非法组合）。
- `peaks`（必填，1～200 个）与 `mass_tolerance`（必填，>0）：含义与校验规则同单标签接口。

拟合模型与单标签接口一致，只是两个被标记元素的分布同时替换为各自的富集分布；除两个被拟合元素外，其余元素的天然包络与两个候选丰度都无关，只计算一次并在整个二维网格中复用；第一个元素的富集多项式在外层丰度固定时同样只计算一次。

响应：

```json
{
  "formula": "C2H5NO2",
  "best_abundances": [
    {"label": "13C", "best_abundance_percent": 30.0},
    {"label": "15N", "best_abundance_percent": 70.0}
  ],
  "mse": 0.0,
  "search_points": 10201,
  "pair_count": 12,
  "pairs": [
    {"observed_mass": 75.032028, "observed_intensity": 42.86, "predicted_mass": 75.032028, "predicted_intensity": 42.86}
  ]
}
```

- `best_abundances`：两个标签各自的最佳丰度（百分比），顺序与请求中的 `labels` 一致。
- `search_points`：二维搜索点数，即两个标签网格点数之积（如各自 1.0 步长时为 101×101 = 10201；两个网格都始终包含 100% 端点）。
- `mse`、`pair_count`、`pairs`：含义与单标签接口相同，配对规则（容差内最近、一对一、未匹配峰按强度 0 计）也完全一致。
- 平局规则：MSE 相同时，取两项丰度按请求顺序（第一标签优先、第二标签其次）升序的字典序最小解，等价于"第一标签丰度尽量低；仍相同则第二标签丰度尽量低"。

### 错误定位

非法分子式或标记返回 HTTP 422，并指出出错字符的 **0 基位置**：

```json
{
  "detail": {
    "error": "invalid_formula",
    "message": "Element 'C' has no known isotope with mass number 20; check the label mass",
    "position": 2
  }
}
```

可检测的错误包括：未知元素符号、元素计数为 0、未闭合的 `[`、标记格式错误、未知标记质量数、孤立数字、非法字符等。

`/api/v1/enrichment-fit` 的 422 错误同样定位明确：

```json
{
  "detail": {
    "error": "invalid_label",
    "message": "Element 'C' has no known isotope with mass number 20; check the label mass",
    "position": 0
  }
}
```

`/api/v1/dual-enrichment-fit` 中两个标签属于同一元素时返回可定位的 422（`index` 为标签列表中出错标签的 0 基索引，指向第二个标签）：

```json
{
  "detail": {
    "error": "invalid_labels",
    "message": "The two labels must target different elements, but '12C' targets element 'C' already named by '13C'",
    "index": 1
  }
}
```

- `invalid_formula`：分子式非法（`position` 为分子式中的 0 基字符位置）。
- `invalid_adduct`：`/api/v1/ion-envelope` 的加合物非法——片段缺少 `+`/`-` 前导符、系数为 0、片段为空/含空白、片段内分子式错误（`position` 为加合物字符串中的 0 基字符位置；片段内错误的位置已换算为绝对位置），或应用片段后某元素/标记计数变为负数（`position` 指向导致首次负计数的片段首字符）：

```json
{
  "detail": {
    "error": "invalid_adduct",
    "message": "Adduct fragment would drive element 'H' count from 2 to -1; counts must not become negative",
    "position": 2
  }
}
```

  注意：`/api/v1/ion-identify` 与 `/api/v1/ion-mixture-fit` 对候选加合物**不**返回上述 422——非法候选以 `valid: false` 出现在 200 响应的 `candidates` 中，并携带其请求 `index` 与同一套 `error`/`message`/`position`；其中 mixture-fit 还会跳过所有涉及该候选的候选对（有效候选少于两个时 `best_mixture` 为 `null`）。仅请求级输入（分子式、峰列表、参数边界）非法时才返回 422。

- `invalid_label`：标签非法——格式错误（应为 `13C` 形式）、未知元素符号、未知质量数（`position` 为标签字符串中的 0 基字符位置）。
- `invalid_labels`：双标签组合非法——两个标签属于同一元素（`index` 为标签列表中出错标签的 0 基索引）。
- `invalid_peaks`：观测峰列表为空、质量非递增或重复（`index` 为出错峰的 0 基索引，空列表时为 `null`）。
- 参数越界（`step_percent` 超出 0.5～10、`mass_tolerance` 非正、`labels` 不是恰好两个、峰数超过 200、强度为负等）由请求校验返回标准 422 并指明字段位置。

### 其他端点

- `GET /api/v1/elements`：列出全部支持的元素及其天然同位素质量/丰度。
- `GET /health`：健康检查。

## 配置（环境变量）

均可在 `docker-compose.yml` 的 `environment` 段覆盖：

| 变量 | 默认值 | 说明 |
|---|---|---|
| `MSPC_DB_PATH` | 容器内 `/data/mspc.db`，本地 `./mspc.db` | SQLite 数据库文件路径 |
| `MSPC_MERGE_TOLERANCE` | `0.0001` | 峰合并质量容差（Da） |
| `MSPC_PRUNE_THRESHOLD` | `1e-12` | 卷积过程中按相对强度剪枝的阈值 |
| `MSPC_TOP_N_MAX` | `200` | `top_n` 上限 |
| `MSPC_TOP_N_DEFAULT` | `20` | `top_n` 默认值 |

## 本地开发（不用 Docker）

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --reload          # 服务: http://127.0.0.1:8000
pytest                                 # 运行测试
```

## 实现说明

- **多项式卷积**：每种元素的同位素模式表示为多项式 `Σ aᵢ·x^mᵢ`（丰度 × 质量的幂），分子包络为各元素多项式（按原子数取幂，使用快速幂）之积；多项式乘法即 (质量, 强度) 峰列表的卷积。纯同位素标记贡献固定质量位移。
- **峰合并**：排序后相邻质量差 < 0.0001 Da 的峰融合为一个，质量取强度加权平均，强度取和。
- **离子包络**：`POST /api/v1/ion-envelope` 将加合物（`app/adduct.py`，形如 `+2H-H2O` 的连续片段，系数乘整个片段）施加到解析后的原子计数上，用 `|charge| × MSPC_MERGE_TOLERANCE` 的中性质量容差卷积，再按 `(mass − charge·mₑ)/|charge|`（`mₑ = 0.000548579909065 Da`）逐峰换算为 m/z。
- **离子识别**：`POST /api/v1/ion-identify`（`app/identify.py`）对 1～20 个候选（加合物, 电荷）逐个复用同一套加合物解析、离子包络计算与基峰归一化，观测峰同样归一化一次后复用；每个候选用与富集拟合相同的 `pair_and_score` 做容差内质量差全局从小到大的一对一贪心匹配，未匹配峰以零强度计入归一化 MSE。有效候选按 MSE 升序、输入索引作稳定平局键排序；非法候选捕获既有 `AdductError`（携带位置）后以 `valid: false` 附在结果末尾，不影响其他候选评分。
- **双离子混合拟合**：`POST /api/v1/ion-mixture-fit`（`app/mixture.py`）对 2～10 个候选（加合物, 电荷）复用离子识别的预测与归一化；枚举每个不同有效候选的无序对，在首候选 0%～100% 的比例网格上，将两个保持各自基峰归一化的包络按比例缩放叠加、用 `mass_tolerance` 经既有 `merge_peaks` 融合（强度加权均值质量、强度求和）后再次基峰归一化，再用同一套 `pair_and_score` 与观测峰一对一配对并把双侧未匹配峰计入 MSE。每对取网格上 MSE 最小的比例（浮点噪声内的平局取较低首候选比例），候选对再按 MSE 升序、较小索引对作平局键排名。非法候选返回其位置错误并跳过涉及它的全部组合；有效候选少于两个时无组合、最佳结果为空。
- **单同位素质量**：按 IUPAC 惯例取各元素天然丰度最高的同位素质量之和。
- **富集拟合**：对每个候选丰度 p（0%～100% 按步长网格搜索），将被标记元素的同位素分布替换为"标记同位素占 p、其余天然同位素按比例分享 1−p"，复用同一套卷积生成归一化包络（基峰=100，观测峰同样归一化）；观测峰与容差内最近的预测峰贪心配对，未匹配峰强度按 0 计，以全部配对的强度均方误差最小者为最佳丰度，同分取较低丰度。除被拟合元素外其余元素的包络与候选丰度无关，只计算一次。
- **双元素联合拟合**：`POST /api/v1/dual-enrichment-fit` 在两个不同元素标签的丰度上做二维网格搜索，组合两种富集分布后沿用同一套卷积、一对一配对与含未匹配峰的 MSE。与两个候选丰度都无关的"其余元素包络"只卷积一次并在全网格复用；第一元素的富集多项式在其外层丰度固定时只计算一次。MSE 相同按两项丰度（请求顺序）升序字典序取解，`search_points` 为两网格点数之积。
- **同位素数据**：`app/isotopes.py` 内置约 50 种常见元素的 AME/CIAAW 同位素质量与天然丰度，可按需扩充。
- **持久化**：SQLite 存计算结果缓存（`envelope_cache`）与请求日志（`request_log`），schema 随服务启动自动创建。
