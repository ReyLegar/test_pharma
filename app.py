import streamlit as st
import numpy as np
from scipy.integrate import odeint

# --- Конфигурация страницы ---
st.set_page_config(page_title="PharmaFunc MVP", layout="wide")
st.title("PharmaFunc Dynamics MVP: Коррекция дозы ванкомицина")

# --- Базовые параметры препарата (JSON-база данных) ---
drug_db = {
    "Vancomycin": {
        "CL_population": 4.5,  # Средний клиренс (л/ч) при СКФ=90 мл/мин
        "Vd_population": 0.7,  # Объем распределения (л/кг)
        "therapeutic_range": (15, 25)  # Терапевтическое окно (мкг/мл)
    }
}

# --- Ввод данных пациента ---
with st.sidebar:
    st.header("Параметры пациента")
    weight = st.number_input("Вес (кг)", min_value=30, max_value=150, value=70)
    age = st.number_input("Возраст", min_value=18, max_value=100, value=45)
    scr = st.number_input("Креатинин сыворотки (мг/дл)", min_value=0.1, max_value=10.0, value=1.0)
    drug = st.selectbox("Препарат", list(drug_db.keys()))

# --- Расчет СКФ по формуле CKD-EPI ---
def calculate_gfr(age, scr, male=True):
    if male:
        kappa = 0.9
        alpha = -0.411
    else:
        kappa = 0.7
        alpha = -0.329
    gfr = 141 * (kappa ** alpha) * (0.993 ** age) * (1.018 if not male else 1) / scr
    return gfr

gfr = calculate_gfr(age, scr)
st.sidebar.write(f"**СКФ (мл/мин/1.73 м²):** {gfr:.1f}")

# --- PK-модель (однокомпартментная) ---
def pk_model(C, t, Dose, CL, Vd, ka, ke):
    dCdt = (ka * Dose) / Vd - ke * C
    return dCdt

# --- Коррекция клиренса по СКФ ---
def adjust_cl_for_gfr(cl_population, gfr):
    return cl_population * (gfr / 90)

# --- Расчет дозы ---
def calculate_dose(drug_name, weight, gfr):
    drug_params = drug_db[drug_name]
    CL_population = drug_params["CL_population"]
    Vd = drug_params["Vd_population"] * weight
    
    # Индивидуальный клиренс
    CL = adjust_cl_for_gfr(CL_population, gfr)
    
    # Целевая концентрация (средняя терапевтическая)
    C_target = np.mean(drug_params["therapeutic_range"])
    
    # Стационарная концентрация: C_ss = (Dose * F) / (CL * tau)
    # Предполагаем F=1 (биодоступность) и tau=24 часа
    tau = 24
    Dose = (C_target * CL * tau) 
    
    return Dose, CL, Vd

# --- Вывод результатов ---
if st.button("Рассчитать дозу"):
    Dose, CL, Vd = calculate_dose(drug, weight, gfr)
    
    st.subheader("Рекомендация")
    st.write(f"**Клиренс препарата:** {CL:.2f} л/ч")
    st.write(f"**Объем распределения:** {Vd:.2f} л")
    st.write(f"**Расчетная доза:** {Dose:.2f} мг каждые 24 часа")
    
    # Визуализация концентрации
    t = np.linspace(0, 24, 100)
    C0 = 0  # Начальная концентрация
    ka = 1.2  # Константа всасывания (предположение)
    ke = CL / Vd  # Константа выведения
    
    solution = odeint(pk_model, C0, t, args=(Dose, CL, Vd, ka, ke))
    st.line_chart(solution, use_container_width=True)
    st.caption("Прогноз концентрации препарата в крови за 24 часа.")

# --- Disclaimer ---
st.markdown("---")
st.caption("⚠️ Это упрощенная демонстрация. Не используйте для реального лечения!")