# SQL

> sample_basic_information

```sql
CREATE TABLE sample_basic_information (
    sample_id VARCHAR(255) PRIMARY KEY,
    date DATE,
    operator VARCHAR(255),
    well_name VARCHAR(255),
    location VARCHAR(255),
    field VARCHAR(255)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;
```

> sample_composition_information

```sql
CREATE TABLE sample_composition_information (
    id INT AUTO_INCREMENT PRIMARY KEY,
    sample_id VARCHAR(255) NULL,
    na_aq DOUBLE DEFAULT 0,
    k_aq DOUBLE DEFAULT 0,
    mg_aq DOUBLE DEFAULT 0,
    ca_aq DOUBLE DEFAULT 0,
    sr_aq DOUBLE DEFAULT 0,
    ba_aq DOUBLE DEFAULT 0,
    feii_aq DOUBLE DEFAULT 0,
    zn_aq DOUBLE DEFAULT 0,
    pb_aq DOUBLE DEFAULT 0,
    cl_aq DOUBLE DEFAULT 0,
    so4_aq DOUBLE DEFAULT 0,
    f_aq DOUBLE DEFAULT 0,
    br_aq DOUBLE DEFAULT 0,
    si_aq DOUBLE DEFAULT 0,
    alk_bicarbonate_aq DOUBLE DEFAULT 0,
    alk_carbonate_aq DOUBLE DEFAULT 0,
    orgacid_acetate_aq DOUBLE DEFAULT 0,
    ammonia_aq DOUBLE DEFAULT 0,
    b_aq DOUBLE DEFAULT 0,
    tds_aq DOUBLE DEFAULT 0,
    density_stp DOUBLE NULL,
    co2_pct_g DOUBLE DEFAULT 0,
    option_use_h2sg INT DEFAULT 0,
    h2s_pct_g DOUBLE DEFAULT 0,
    h2s_aq DOUBLE DEFAULT 0,
    ph_stp DOUBLE NULL,
    q_gas DOUBLE DEFAULT 0,
    q_oil DOUBLE DEFAULT 0,
    q_water DOUBLE DEFAULT 1000
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

```

> temperature_pressure_condition_information

```sql
CREATE TABLE temperature_pressure_condition_information (
    id INT AUTO_INCREMENT PRIMARY KEY,
    sample_id VARCHAR(255) NULL,
    t_initial DOUBLE DEFAULT 77,
    t_final DOUBLE DEFAULT 77,
    p_initial DOUBLE DEFAULT 14.7,
    p_final DOUBLE DEFAULT 14.7,
    api DOUBLE NULL,
    sg_g DOUBLE NULL,
    q_meoh DOUBLE DEFAULT 0,
    q_meg DOUBLE DEFAULT 0,
    strongacid_aq DOUBLE DEFAULT 0,
    strongbase_aq DOUBLE DEFAULT 0,
    conc_multiplier DOUBLE NOT NULL,
    t_ph DOUBLE NULL,
    p_ph DOUBLE NULL,
    t_q DOUBLE NULL,
    p_q DOUBLE NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

```

> sample_oil_phase_information

```sql
CREATE TABLE sample_oil_phase_information (
    id INT AUTO_INCREMENT PRIMARY KEY,
    sample_id VARCHAR(255) NULL,
    c1_o DOUBLE NULL,
    co2_o DOUBLE NULL,
    h2s_o DOUBLE NULL,
    c2_o DOUBLE NULL,
    c3_o DOUBLE NULL,
    ic4_o DOUBLE NULL,
    nc4_o DOUBLE NULL,
    ic5_o DOUBLE NULL,
    nc5_o DOUBLE NULL,
    c6_o DOUBLE NULL,
    c7_c12_o DOUBLE NULL,
    c13_c25_o DOUBLE NULL,
    c26_c80_o DOUBLE NULL,
    n2_o DOUBLE NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

```

>configure_parameter_information

```sql
CREATE TABLE configure_parameter_information (
    sample_id VARCHAR(255) NULL,
    Option_Alk INT DEFAULT 0,
    Option_Defined_TP INT DEFAULT 0,
    Option_TP_for_pH INT DEFAULT 0,
    Option_TP_for_Q INT DEFAULT 0,
    Option_EoS INT DEFAULT 0,
    Option_Water_HC INT DEFAULT 0
);
```

>newssp_software_data_dictionary

```sql
CREATE TABLE newssp_software_data_dictionary (
    Field_Name VARCHAR(200) NOT NULL,       
    English_Name VARCHAR(200),              
    Chinese_Name VARCHAR(200),              
    Unit VARCHAR(200),                       
    PRIMARY KEY (Field_Name)
);


INSERT INTO newssp_software_data_dictionary (Field_Name, English_Name, Chinese_Name, Unit) VALUES
('Date', '', '采样日期', ''),
('Operator', '', '采样人员', ''),
('Well Name', '', '油井名称', ''),
('Location', '', '所在位置', ''),
('Field', '', '油田名称', ''),
('Option_Alk', 'Four calc options', '碳酸体系计算选项', '0-CO2%+Alk,1-pH+Alk,2-CO2%+pH,3-HCO3+CO3'),
('Option_Defined_TP', 'Use TP on Calcite sheet?', '是否使用calcite页面的温压值', '1-Yes / 0-No'),
('Option_TP_for_pH', 'T, P for pH', '是否使用输入的pH测量时的温压值', '1-T,P;0-STP'),
('Option_TP_for_Q', 'T,P for G/O/W', '是否使用输入的流量测量时的温压值', '1-T,P;0-STP'),
('Option_EoS', 'Use Flash Calculator', '是否使用flash calculator', '0-No;1-Res Comp;2-Res Comp+pH or PCO2;3-Psuedo Comp+pH or PCO2'),
('Option_Water_HC', 'Water saturation in HC phase?', '是否假设水在有机相中饱和', '1-No H2O/0-With H2O'),
('Na_aq', 'Na+', '', '(m)'),
('K_aq', 'K+', '', '(m)'),
('Mg_aq', 'Mg2+', '', '(m)'),
('Ca_aq', 'Ca2+', '', '(m)'),
('Sr_aq', 'Sr2+', '', '(m)'),
('Ba_aq', 'Ba2+', '', '(m)'),
('FeII_aq', 'Fe2+', '', '(m)'),
('Zn_aq', 'Zn2+', '', '(m)'),
('Pb_aq', 'Pb2+', '', '(m)'),
('Cl_aq', 'Cl-', '', '(m)'),
('SO4_aq', 'SO42-', '', '(m)'),
('F_aq', 'F-', '', '(m)'),
('Br_aq', 'Br-', '', '(m)'),
('Si_aq', 'Silica', '', '(m as Si)'),
('Alk_Bicarbonate_aq', 'Total Alkalinity', '总碱度', '(m as HCO3)'),
('Alk_Carbonate_aq', 'CO32- Alkalinity', 'CO32-碱度', '(m as CO3)'),
('OrgAcid_Acetate_aq', 'Carboxylates', '有机酸', '(m as Acetate)'),
('Ammonia_aq', 'Ammonia', '', '(m as NH3)'),
('B_aq', 'Borate', '', '(m as B)'),
('TDS_aq', 'TDS (Measured)', '测量的TDS', '(mg/l)'),
('Density_STP', 'Calc. Density (STP)', 'STP条件下密度', '(g/ml)'),
('CO2_pct_g', 'CO2 Gas Analysis', 'CO2气相浓度', '(%)'),
('Option_Use_H2Sg', 'Use H2S Gas Analysis', '使用H2S气相浓度（1）或H2S溶解浓度（0）', '1-Yes / 0-No'),
('pH_STP', 'pH, measured', 'pH在STP下的测量值', 'pH'),
('Q_Gas', 'Gas/day(thousand cf/day)', '气体流量', '(Mcf/D)'),
('Q_Oil', 'Oil/Day', '油流量', '(B/D)'),
('Q_Water', 'Water/Day', '水流量', '(B/D)'),
('T_initial', 'Initial T', '起点温度', '(°F)'),
('T_final', 'Final T', '终点温度', '(°F)'),
('P_initial', 'Initial P', '起点压力', '(psia)'),
('P_final', 'Final P', '终点压力', '(psia)'),
('API', 'API Oil Grav.', '油API密度', 'API grav.'),
('SG_g', 'Gas Sp.Grav.', '气体比重', 'Sp.Grav.'),
('Q_MeOH', 'MeOH/Day', 'MeOH流量', '(B/D)'),
('Q_MEG', 'MEG/Day', 'MEG流量', '(B/D)'),
('StrongAcid_aq', 'H+ (Strong acid)', '加入强酸量', '(m)'),
('StrongBase_aq', 'OH- (Strong base)', '加入强碱量', '(m)'),
('Conc_Multiplier', 'Conc. multiplier', '浓度系数', ''),
('T_pH', 'Temp. for pH meas.', 'pH测量温度', '(°F)'),
('P_pH', 'Pres. for pH meas.', 'pH测量压力', '(psia)'),
('T_Q', 'T for fluids meas.', '流量测试温度', '(°F)'),
('P_Q', 'P for fluids meas.', '流量测试压力', '(psia)'),
('C1_o', 'C1', '', 'mole %'),
('CO2_o', 'CO2', '', 'mole %'),
('H2S_o', 'H2S', '', 'mole %'),
('C2_o', 'C2', '', 'mole %'),
('C3_o', 'C3', '', 'mole %'),
('iC4_o', 'iC4', '', 'mole %'),
('nC4_o', 'nC4', '', 'mole %'),
('iC5_o', 'iC5', '', 'mole %'),
('nC5_o', 'nC5', '', 'mole %'),
('C6_o', 'C6', '', 'mole %'),
('C7_C12_o', 'C7-C12', '', 'mole %'),
('C13_C25_o', 'C13-C25', '', 'mole %'),
('C26_C80_o', 'C26-C80', '', 'mole %'),
('N2_o', 'N2', '', 'mole %');


```
