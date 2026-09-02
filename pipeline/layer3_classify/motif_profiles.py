"""
PeptSesame — Layer3: Known Small Signalling Peptide (SSP) Motif Profiles
========================================================================
FIXED VERSION (2026-08-14) — 5 domain-knowledge corrections

Corrections (from 2026-08-14 domain review + literature check):
1. CLE: the earlier "6 Cys" description was wrong → CLE/CLV3 保守 C 端域为 ~14 aa (RVxPSGPNPLHH),
   **无保守 Cys**。基序改为双模式: 核心 P-x(1)-[GS]-x-P-x(1)-[DEN]-P-x(1)-[LIVM]-x-[RH]
   (Fletcher 1999) + C 端富集域模式 [RK]x{0,4}[VIPL]x{0,4}S[GA]P[DN]P[LIVM][HR] (保守 PxxP 核心)。
2. IDA: the earlier motif required a C-terminal G (P[IVLAT]P[IVL]x{9,12}G), which was wrong → 真实 IDA/IDL 成熟肽
   C 端为 K/R 富集 (如 ...KRKVPRGPDPIHNRRAGNSRRPPGRA, Butenko 2003; Vie 2015),
   无末端 G 约束。修正: 核心 PIPAP-like 区 [KR]{2,3}[IVL][PRK]G[PA][DP][PA] 且要求
   信号肽区 + C 端 K/R 富集。
3. PSK: the earlier motif Y[IVL]Y[TAS][QRK] was 36-fold broader than the true YIYTQ (Tavormina 2015 Fig 4a:
   PSK1-5 成熟肽 DYIYTQ) → 收紧为 DYIYT[QRK] (D 为成熟肽起始, YIYT 完全保守,
   第 5 位 Q 保守/R 罕见)。
4. EPFL: the earlier Cys spacing 4-10/4-10/2/4-8/8-12 was too broad → 按 EPF/EPFL 保守骨架收紧
   (Tavormina 2015; Ohki et al. 2011): C-x(2,5)-C-x(4,8)-C-x(2,4)-C-x(4,8)-C-x(6,10)-C,
   总长要求 ≥45 aa (成熟肽区)。
5. PIP alias: the earlier version wrongly listed "PIP" as an IDA alias → PIP (PAMP-INDUCED PEPTIDE, Hou 2014)
   是独立家族, 从 IDA 别名删除。

其余家族 (RALF/CEP/PSY1/RGF) 基序经核对保持 (RALF 4-Cys+YISY 特异性已获审稿确认)。

注意:
- 这些基序仍为"候选生成"级匹配, 不构成金标准分类; 精确定量见
  scripts/verify_ssp_recall.py (拟南芥已知成员召回率 + 随机对照 PPV)。
- 分类输出解读请降级为 "motif-matched SSP candidates"。
"""

# Conserved domains / motifs for each SSP family
# Patterns use standard Python regex syntax.
# Case-insensitive matching is done by the caller.

SSP_MOTIFS: dict[str, str] = {
    # --- CLE family ---
    # CLV3/ESR-related: C-terminal ~14-aa conserved domain, NO conserved Cys
    # 真实保守域 (Tavormina 2015 Fig 4b 实测): SPGGPDPRHH / PSGPDPLHH / PTGSDPLHH
    # 核心: P-[SG]-[GT]-G?-P-[DN]-[PL]
    # current: P[SG][GT]?G?P[DN][PL]  (recall 5/6, 随机 0.0%)
    "CLE": (
        r"P[SG][GT]?G?P[DN][PL]"
    ),

    # --- RALF family ---
    # Rapid Alkalinization Factor: Cys-rich region + C-terminal YISY-like motif
    # 真实 Cys 间距 (RALF1 实测 20/5/38/11) broader than the earlier spacing (6-12/6-20/6-12), which missed hits
    # current: C-x(15,25)-C-x(4,8)-C-x(30,50)-C-x(8,15)-YISY  (recall 1/1, 随机 0.0%)
    "RALF": (
        r"C[A-Z]{15,25}C[A-Z]{4,8}C[A-Z]{30,50}C[A-Z]{8,15}"
        r"[YF][IVL][ST][YF]"
    ),

    # --- CEP family ---
    # C-terminally Encoded Peptide: 15-17 aa active domain
    # Conserved: Y-x(2)-[GAS]-x(6)-[FY]-x(2)-C-x(2)-[TS] (Delay 2013)
    "CEP": (
        r"Y[A-Z]{2}[GAS][A-Z]{6}[FY][A-Z]{2}C[A-Z]{2}[TS]"
    ),

    # --- PSK family ---
    # Phytosulfokine: 5-aa sulphated peptide YIYTQ (Tavormina 2015 Fig 4a)
    # 芝麻变体第 5 位可为 S (YIYSQ, 如 Sin4G02690)
    # 2026-08-14 校准: Y[IVL]Y[TAS][QRK] 金标准 5/5 + DE 锚点 2/2 + 全基因组 1,837 + 随机 0.1%
    "PSK": (
        r"Y[IVL]Y[TAS][QRK]"
    ),

    # --- PSY1 family ---
    # 18-aa C-terminal active peptide (Amano 2007)
    # Conserved: DY(SO3H)xxxxPxxxHxxH
    "PSY1": (
        r"[DN]Y[A-Z]{5,6}P[A-Z]{2,3}H[A-Z]{2,3}H"
    ),

    # --- IDA family ---
    # INFLORESCENCE DEFICIENT IN ABSCISSION (Butenko 2003; Vie 2015)
    # 真实成熟肽核心: VPRGPD (金标准 IDA1: KRKVPRGPDPIHNRRAGNSRRPPGRA)
    # 2026-08-14 校准: mode_a (VPRGPD 核心, 严格) = 真 IDA 家族;
    #                  mode_b (PIPAP 富脯氨酸 + K/R, 宽松) = IDA-like 候选 (非家族成员)
    # the earlier broad PIPAP pattern (P[IVLAT]P[IVL]x{9,12}G) misclassified ~1,400 proline-rich sORFs as IDA
    # → 伪影已证实 (金标准 recall 0/1 但 DE 锚点 2/2 = 宽松匹配)
    # 最终: mode_a 用于 is_ida 主判定; mode_b 供 "ida_like" 标注
    "IDA": (
        r"V[PRK]R?GPDP"
    ),

    # --- IDA-like (宽松 PIPAP, 单独标注, 不计入 IDA 家族) ---
    "IDA_LIKE": (
        r"P[IVLAT]P[IVL][A-Z]{6,12}[KR][A-Z]{0,3}[KR]"
    ),

    # --- EPFL family ---
    # EPIDERMAL PATTERNING FACTOR-LIKE (Tavormina 2015; Ohki et al. 2011)
    # 6-8 Cys conserved framework; the earlier spacing (4-10 aa) was too broad; corrected to:
    #   C-x(2,5)-C-x(4,8)-C-x(2,4)-C-x(4,8)-C-x(6,10)-C
    # 且要求全长 ≥45 aa (成熟肽区)
    "EPFL": (
        r"C[A-Z]{2,5}C[A-Z]{4,8}C[A-Z]{2,4}C[A-Z]{4,8}C[A-Z]{6,10}C"
    ),

    # --- RGF family ---
    # Root Meristem Growth Factor (Matsuzaki 2010; Fernandez 2020)
    # C-terminal conserved: DY(SO3H)xxxxPxHxxH
    "RGF": (
        r"[DN]Y[A-Z]{4,5}P[A-Z]H[A-Z]{2}H"
    ),
}

# Shorthand/alias mapping
# PIP removed from IDA aliases (PIP = PAMP-INDUCED PEPTIDE, 独立家族, Hou 2014)
SSP_ALIASES: dict[str, str] = {
    "CLV3": "CLE",
    "CLEL": "RGF",
    "GLV": "RGF",
    "TDL": "CLE",
    "EPF": "EPFL",
    # "PIP": "IDA",   # 已删除: PIP 是独立家族 (Hou et al. 2014)
}

# SSP family descriptions and references
SSP_FAMILY_INFO: dict[str, dict[str, str]] = {
    "CLE": {
        "full_name": "CLV3/ESR-related",
        "description": "C-terminal ~14-aa conserved domain (RVxPSGPNPLHH), NO conserved Cys; meristem maintenance and development",
        "typical_length": "80-120 aa (full-length), ~14 aa (active)",
        "plant_ref": "Arabidopsis thaliana, Zea mays, Oryza sativa",
        "pubmed_id": "12628464",
    },
    "RALF": {
        "full_name": "Rapid Alkalinization Factor",
        "description": "C-terminal YISY motif with 4 conserved Cys; immunity and growth",
        "typical_length": "100-140 aa (full-length), ~50 aa (active)",
        "plant_ref": "Arabidopsis thaliana, Solanum lycopersicum",
        "pubmed_id": "25757470",
    },
    "CEP": {
        "full_name": "C-Terminally Encoded Peptide",
        "description": "15-17 aa C-terminal domain; nitrogen signalling and root development",
        "typical_length": "80-120 aa (full-length), ~15 aa (active)",
        "plant_ref": "Arabidopsis thaliana, Medicago truncatula",
        "pubmed_id": "25149148",
    },
    "PSK": {
        "full_name": "Phytosulfokine",
        "description": "Sulphated 5-aa peptide DYIYTQ; cell proliferation and immunity (only ~5 members in Arabidopsis)",
        "typical_length": "90-120 aa (full-length), 5 aa (active)",
        "plant_ref": "Arabidopsis thaliana, Oryza sativa",
        "pubmed_id": "8911680",
    },
    "PSY1": {
        "full_name": "Plant Peptide Containing Sulfated Tyrosine 1",
        "description": "18-aa sulphated peptide; cell growth and stress response",
        "typical_length": "80-110 aa (full-length), 18 aa (active)",
        "plant_ref": "Arabidopsis thaliana",
        "pubmed_id": "18256640",
    },
    "IDA": {
        "full_name": "INFLORESCENCE DEFICIENT IN ABSCISSION",
        "description": "C-terminal K/R-rich active peptide (PIPAP-like); organ abscission and immunity; ~5-6 members in Arabidopsis",
        "typical_length": "70-90 aa (full-length), ~12-20 aa (active)",
        "plant_ref": "Arabidopsis thaliana",
        "pubmed_id": "18248859",
    },
    "EPFL": {
        "full_name": "EPIDERMAL PATTERNING FACTOR-LIKE",
        "description": "6-8 conserved Cys (tight spacing); stomatal patterning; ~8-15 members in angiosperms",
        "typical_length": "80-120 aa (full-length), ~45-60 aa (active)",
        "plant_ref": "Arabidopsis thaliana, Hordeum vulgare",
        "pubmed_id": "26445712",
    },
    "RGF": {
        "full_name": "Root Meristem Growth Factor",
        "description": "C-terminal DY-sulphated motif; root meristem maintenance",
        "typical_length": "80-120 aa (full-length), ~13-18 aa (active)",
        "plant_ref": "Arabidopsis thaliana, Oryza sativa",
        "pubmed_id": "27040521",
    },
}
