# Field-bias volume: UCD_6 high-confinement boundary-value audit
# 場修正在 UCD_6 高圍壓案例的邊值試驗

Companion to `FIELD_BIAS_VOLUME_VALIDATION.md`. That note states the
high-confinement window (admitted mean pressure above the 40-kPa transition)
is **not calibrated**. This is a measured result inside exactly that window.
本檔為 `FIELD_BIAS_VOLUME_VALIDATION.md` 的補充。該文件註明「高於 40 kPa
轉換壓力的高圍壓窗尚未校正」；本檔即該窗內的一次實測。

First: thank you for adopting the serialized reversal latch and `-noBiasVolume`.
Both goldens pass here (pAnchor `106.11646292406457803281` bit-identical,
revision-4 save/restore with serialized latch/field-bias continuity).
先致謝：您已納入序列化 reversal latch 與 `-noBiasVolume`；此處兩個 golden 均通過
（pAnchor 逐位元相同、revision-4 含 latch/field-bias 序列化連續性）。

---

## Test / 試驗

LEAP-Asia-2019 **UCD_6** centrifuge slice, 125×25 SSPbrickUP, one frozen
Ottawa F-65 row, three-stage workflow (`slice_gravity_wenyang.tcl`), reversal
latch on. **Only one thing changed** between the two runs, exactly as
suggested: `-noBiasVolume` → `-fieldBiasVolume`.
LEAP-Asia-2019 **UCD_6** 離心機 slice，125×25 SSPbrickUP，同一組凍結 Ottawa F-65，
三階段工作流、latch 開。兩次分析**只差一處**：`-noBiasVolume` 改成 `-fieldBiasVolume`。

UCD_6 is dense (Dr 0.666); its pore-pressure stations sit at σv0′ ≈ 39 / 30 /
20 / 9 kPa near surface but the deep column reaches **178 kPa**, well above the
40-kPa transition — so `-fieldBiasVolume` is **active** here, unlike RPI-A
(max admitted mean ≈ 28 kPa, correction inactive).
UCD_6 為密砂（Dr 0.666）；深層柱達 **178 kPa**，遠高於 40 kPa 轉換壓力，故
`-fieldBiasVolume` 在此**作用**；而 RPI-A 最大容許平均壓約 28 kPa，修正不作用。

---

## 1. Analysis time — the headline / 分析時間（最關鍵）

Same mesh, same nSub (8), same everything but the one flag:
同網格、同 nSub（8）、除該旗標外全同：

| run | flag | wall time | recovered failures |
|---|---|---:|---:|
| W9N | `-noBiasVolume` (mechanism OFF) | **62 min** | **0** |
| W10 | `-fieldBiasVolume` (mechanism ON + bounded field correction) | **240 min** | **131** |

**≈ 3.9× slower.** The cost is the live mechanism: with
`bias_reversible_volume` active, every genuine reversal near zero effective
stress (the liquefied loose/medium elements) fires the volume term, the coupled
Newton step struggles, dt is cut and re-tried — 131 recovered cuts vs 0. The
latch keeps it from diverging (no DNF), but the near-liquefaction grind returns.
**約慢 3.9 倍。** 成本來自機制活著：`bias_reversible_volume` 一活，近零有效應力
（液化中的鬆／中密元素）每次真實反轉都觸發體積項，耦合 Newton 收斂困難、砍 dt 重試
——131 次回復失敗 vs 0。latch 讓它不至於 DNF，但近液化的 grind 回來了。

This matches the wider pattern we measured: ungated (mechanism on, no latch) =
5733 failures / DNF; `-noBiasVolume` = 0–16 failures / fastest; field-corrected
+ latch = 131 failures / 3.9× wall.
這與我們量到的整體規律一致：無閘（機制開、無 latch）＝5733 失敗／DNF；
`-noBiasVolume`＝0–16 失敗／最快；場修正＋latch＝131 失敗／3.9× wall。

---

## 2. Pore pressure — overshoot then over-dissipation / 水壓：先過衝後過度消散

r_u,max and retained r_u at t = 25 s vs the centrifuge:
r_u 峰值與 t=25 s 的殘餘值，對照離心機：

| station | σv0′ (kPa) | r_u,max cen / W9N / W10 | r_u@25s cen / W9N / W10 |
|---|---:|---|---|
| P1 | 38.6 | 0.83 / 0.91 / **1.04** | 0.49 / 0.55 / **0.30** |
| P2 | 30.5 | 0.99 / 0.78 / **1.53** | 0.64 / 0.62 / **0.36** |
| P3 | 20.4 | 1.06 / 0.91 / **1.86** | 0.85 / 0.61 / **0.34** |
| P4 | 8.8 | 1.05 / 0.80 / **2.08** | 1.02 / 0.63 / **0.34** |

`-fieldBiasVolume` (W10): during shaking the excess pore pressure oscillates
**well above σv0′** (r_u peaks 1.5–2.1 = large near-zero-σ integration jitter,
the 131 failures), then **over-dissipates** — retained r_u falls to 0.30–0.36,
below both the centrifuge (0.49–1.02) and W9N (0.55–0.63). W9N is smooth and
tracks the centrifuge trend.
`-fieldBiasVolume`（W10）：震動中超額水壓大幅**超過 σv0′**（峰值 1.5–2.1＝近零
有效應力的積分抖動，即那 131 次失敗），震後又**過度消散**——殘餘 r_u 掉到
0.30–0.36，比離心機（0.49–1.02）與 W9N（0.55–0.63）都低。W9N 平滑且貼近離心機趨勢。

---

## 3. Surface displacement / 坡面位移

| run | ux (2 central markers) | ratio | |ln ratio| |
|---|---:|---:|---:|
| centrifuge | 77.5 mm | 1.00 | — |
| W9N `-noBiasVolume` | 71.9 mm | ×0.93 | 0.073 |
| W10 `-fieldBiasVolume` | 86.7 mm | ×1.12 | 0.113 |

W9N is slightly closer to the measurement on UCD_6.
在 UCD_6 上 W9N 略接近量測。

---

## Reading / 解讀

This is **not** a claim that `-fieldBiasVolume` is wrong — it is the measured
behaviour in the high-confinement window your own note flags as uncalibrated.
On UCD_6 (dense, deep σv0′ 178 kPa) the active correction (a) reintroduces the
near-liquefaction convergence cost (3.9× wall, 131 failures), (b) overshoots
peak r_u then over-dissipates, and (c) slightly over-predicts displacement,
whereas `-noBiasVolume` is cleaner and closer.
這**不是**說 `-fieldBiasVolume` 錯了——這是您文件指明「未校正的高圍壓窗」內的實測。
在 UCD_6（密砂、深層 σv0′ 178 kPa），啟動的修正 (a) 帶回近液化收斂成本（3.9× wall、
131 失敗）、(b) r_u 峰值過衝後過度消散、(c) 位移略偏高；而 `-noBiasVolume` 較乾淨、較準。

Suggested next step: exercise `-fieldBiasVolume` on **RPI-A / RPI-B** (low
confinement — correction inactive or weak, the window it was designed for)
before judging it overall; and, if the high-confinement window is to be
supported, the shaking-phase r_u overshoot near zero effective stress looks
like the first thing to bound.
建議下一步：先在**低圍壓的 RPI-A / RPI-B**（修正不作用或很弱，正是它設計服務的窗）
測 `-fieldBiasVolume`，再對它整體下判斷；若要支援高圍壓窗，近零有效應力的震動段
r_u 過衝似乎是首先要綁住的。

Full figures (acc / spectra / pore pressure / displacement, centrifuge vs W9N
vs W10) available on request.
完整圖組（加速度／反應譜／水壓／位移，離心機 vs W9N vs W10）可另附。
