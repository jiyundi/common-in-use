import pandas as pd

obs_manual_adjust = 13 # avg. observed an offset to real payments

# ===== 1. 读取 CSV 文件 =====
file_path = "202505-06.csv"
df = pd.read_csv(file_path, skiprows=2)
df = df[['DATE', 'START TIME', 'USAGE']]
df['datetime'] = pd.to_datetime(df['DATE'] + ' ' + df['START TIME'], 
                                format='%m/%d/%Y %I:%M %p')
df['USAGE_kWh'] = df['USAGE'].astype(float)
df = df[['datetime', 'USAGE_kWh']].sort_values('datetime')

# ===== 2. 基本时间/季节特征 =====
df['month'] = df['datetime'].dt.month
df['hour'] = df['datetime'].dt.hour
df['weekday'] = df['datetime'].dt.weekday  # Monday=0
df['is_summer'] = df['month'].between(5, 9)
df['is_winter'] = ~df['is_summer']

# 峰时段定义
def is_on_peak(row):
    if row['is_summer']:
        return (row['weekday'] < 5) and (15 <= row['hour'] < 19)  # 3–7pm
    else:
        return (row['weekday'] < 5) and ((6 <= row['hour'] < 9) or (18 <= row['hour'] < 21))
df['is_peak'] = df.apply(is_on_peak, axis=1)

# ===== 3. 总用电量 =====
total_kwh = df['USAGE_kWh'].sum()
peak_kwh = df.loc[df['is_peak'], 'USAGE_kWh'].sum()
off_kwh = total_kwh - peak_kwh
peak_kw = df['USAGE_kWh'].max()  # 最大单小时功率视作 kW

# ===== 4. 计算各方案 =====

### (1) Basic Plan ###
def basic_charge(kwh, single_phase=True):
    base_fee = 15 if single_phase else 20
    if kwh <= 500:
        energy = kwh * 0.1304
    elif kwh <= 1000:
        energy = 500 * 0.1304 + (kwh - 500) * 0.1501
    else:
        energy = 500 * 0.1304 + 500 * 0.1501 + (kwh - 1000) * 0.1573
    return base_fee + energy

### (2) Time-of-Use (TOU) ###
def tou_charge(df, single_phase=True):
    base_fee = 12
    total_cost = 0
    for is_summer, group in df.groupby('is_summer'):
        kwh = group['USAGE_kWh'].sum()
        peak = group.loc[group['is_peak'], 'USAGE_kWh'].sum()
        off = kwh - peak
        if is_summer:
            # 夏季档价
            def tier_rate(k):
                if k <= 500:
                    return (peak * 0.1927 + off * 0.1219)
                elif k <= 1000:
                    return (peak * 0.2056 + off * 0.1347)
                else:
                    return (peak * 0.2128 + off * 0.1419)
        else:
            # 冬季档价
            def tier_rate(k):
                if k <= 500:
                    return (peak * 0.1333 + off * 0.1248)
                elif k <= 1000:
                    return (peak * 0.1461 + off * 0.1377)
                else:
                    return (peak * 0.1533 + off * 0.1449)
        total_cost += tier_rate(kwh)
    return base_fee + total_cost

### (3) Peak Demand ###
def peak_demand_charge(df, single_phase=True):
    base_fee = 12 if single_phase else 17
    summer = df.loc[df['is_summer'], 'USAGE_kWh'].sum()
    winter = df.loc[df['is_winter'], 'USAGE_kWh'].sum()
    energy_cost = summer * 0.0883 + winter * 0.0842
    demand_kw = df['USAGE_kWh'].max()
    demand_cost = demand_kw * (11.6 if demand_kw <= 7 else 16.85)
    return base_fee + energy_cost + demand_cost

### (4) Demand TOU ###
def demand_tou_charge(df, single_phase=True):
    base_fee = 12 if single_phase else 17
    total_cost = 0
    for is_summer, group in df.groupby('is_summer'):
        peak = group.loc[group['is_peak'], 'USAGE_kWh'].sum()
        off = group.loc[~group['is_peak'], 'USAGE_kWh'].sum()
        if is_summer:
            energy_cost = peak * 0.1437 + off * 0.0729
        else:
            energy_cost = peak * 0.0843 + off * 0.0758
        total_cost += energy_cost
    demand_kw = df['USAGE_kWh'].max()
    demand_cost = demand_kw * (11.6 if demand_kw <= 7 else 16.85)
    return base_fee + total_cost + demand_cost

# ===== 5. 汇总结果 =====
costs = {
    "Basic":       basic_charge(total_kwh) + obs_manual_adjust,
    "TOU":         tou_charge(df)          + obs_manual_adjust,
    "Peak Demand": peak_demand_charge(df)  + obs_manual_adjust,
    "Demand TOU":  demand_tou_charge(df)   + obs_manual_adjust
}

summary = pd.DataFrame.from_dict(costs, orient='index', columns=['Monthly Cost ($)'])
summary['Monthly Cost ($)'] = summary['Monthly Cost ($)'].round(2)
summary = summary.sort_values('Monthly Cost ($)')

print(summary)
print("\n✅ 最划算方案:", summary.index[0])
