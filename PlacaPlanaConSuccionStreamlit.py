import streamlit as st
from numpy import linalg, zeros, linspace, dot, array, exp
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from matplotlib.pyplot import subplots

# Configuración de la página
st.set_page_config(page_title="Perfil de Velocidades en una Placa Plana con Succión", 
                   page_icon="https://share.google/K4OvrbntlsuX37Y8H.png",
                   layout="wide",
                   initial_sidebar_state="expanded")

def calcular_perfil_velocidades_df(uinf, vw, nu, n, h, esquema):
    """
    Calcula el perfil de velocidades de la placa plana con succión utilizando 
    el esquema de diferencias finitas seleccionado por el usuario. 
    Devuelve el perfil de velocidades resultante.

    Parámetros de entrada:
    - uinf: Velocidad de la corriente libre (float)
    - vw: Velocidad de succión (float)
    - nu: Viscosidad cinemática (float)
    - n: Número de nodos (int)
    - h: Paso espacial (float)
    - esquema: Esquema de diferencias finitas seleccionado (str)

    Parámetro de salida:
    - udf: Perfil de velocidades calculado por diferencias finitas (array)
    """

    # Inicialización de matrices de tamaño n x n en minúscula
    matriz_a = zeros((n, n))
    vector_b = zeros(n)

    # Condiciones de contorno en los extremos de la matriz y el vector
    matriz_a[0, 0] = 1.0
    matriz_a[-1, -1] = 1.0
    vector_b[-1] = uinf

    # Definición de coeficientes según el esquema numérico elegido
    if esquema == "Hacia Atras de 1er orden":
        a = (nu / h**2) - (vw / h)
        b = (-2 * nu / h**2) + (vw / h)
        c = nu / h**2
    elif esquema == "Hacia Adelante de 1er orden":
        a = nu / h**2
        b = (-2 * nu / h**2) - (vw / h)
        c = (nu / h**2) + (vw / h)
    elif esquema == "Centrado de 2do orden":
        a = (nu / h**2) - (vw / (2 * h))
        b = -2 * nu / h**2
        c = (nu / h**2) + (vw / (2 * h))
    else:
        st.error("Esquema de diferencias finitas no implementado.")
        return None

    # Recorrido con un bucle for solo por los nodos internos (de 1 a n-2)
    for i in range(1, n - 1):
        matriz_a[i, i - 1] = a
        matriz_a[i, i]     = b
        matriz_a[i, i + 1] = c

    # Resolución directa del sistema lineal completo
    udf = linalg.solve(matriz_a, vector_b)

    return udf

def objetivo_tiro(a0, ec_diferencial, y, uinf):
    """
    Función objetivo para el método del tiro. Evalúa la diferencia entre la velocidad obtenida 
    al final de la integración y la velocidad de la corriente libre uinf.
    
    Parámetros de entrada:
    - a0: Pendiente inicial (array de un solo elemento)
    - ec_diferencial: Función que define el sistema de EDOs
    - y: Coordenada espacial y (array)
    - uinf: Velocidad de la corriente libre (float)

    Parámetro de salida:
    - error: Diferencia entre la velocidad final y uinf (lista con un elemento)
    """
    # Condición inicial de prueba: u(0) = 0, du/dy(0) = z
    u0 = [0.0, a0[0]]

    # Integramos usando solve_ivp (con RK45) a lo largo de todo el vector y
    sol = solve_ivp(
        ec_diferencial,
        [y[0], y[-1]],
        u0,
        t_eval=y,
        method="RK45",
    )

    # El error es la diferencia entre la velocidad obtenida al final y uinf
    u_final = sol.y[0, -1]
    error = u_final - uinf
    return [error]

def calcular_perfil_velocidades_sm(uinf, vw, nu, y, a_inicial):
    """
    Calcula el perfil de velocidades utilizando el Método del Tiro (Shooting Method).

    Parámetros de entrada:
    - uinf: Velocidad de la corriente libre (float)
    - vw: Velocidad de succión (float)
    - nu: Viscosidad cinemática (float)
    - y: Coordenada y (array de numpy)
    - a_inicial: Valor inicial estimado para la pendiente (float)

    Parámetro de salida:
    - usm: Perfil de velocidades obtenido por el método del tiro (array)
    """
    # Se define la ecuación diferencial de segundo orden como un sistema de primer orden
    ec_diferencial = lambda y_val, u: [u[1], (-vw / nu) * u[1]]

    # fsolve busca la pendiente exacta que hace que el error en el borde superior sea cero
    z_optimo, = fsolve(objetivo_tiro([a_inicial], ec_diferencial, y, uinf), [a_inicial])

    # Con la aceleración óptima ya encontrada, se realiza la integración final
    u0_optimo = [0.0, z_optimo]
    sol_final = solve_ivp(
        ec_diferencial,
        [y[0], y[-1]],
        u0_optimo,
        t_eval=y,
        method="RK45",
    )
    usm = sol_final.y[0]

    return usm

def calcular_perfil_velocidades_analitico(uinf, vw, nu, y):
    """
    Calcula el perfil de velocidades analítico para la placa plana con succión.

    Parámetros de entrada:
    - uinf: Velocidad de la corriente libre (float)
    - vw: Velocidad de succión (float)
    - nu: Viscosidad cinemática (float)
    - y: Coordenada y (array de numpy)

    Parámetro de salida:
    - uan: Perfil de velocidades analítico (array)
    """
    uan = uinf * (1.0 - exp(- (vw / nu) * y))

    return uan

def main():
    """
    Función principal que configura la interfaz de usuario en Streamlit, 
    captura los parámetros de entrada del flujo y numéricos, y coordina 
    la ejecución de los cálculos y gráficos.
    """
    # Presentación de la página y títulos principales
    st.title("Perfil de Velocidades en una Placa Plana con Succión")
    st.subheader("Cátedra de Programación y Análisis Numérico - UTN FRH")
    st.markdown("En esta aplicación se puede obtener el perfil de velocidades sobre una placa plana con succión,")
    st.markdown("que viene determinado por la siguiente ecuación diferencial:")         
    st.latex(r"V_w \frac{du}{dy} + \nu \frac{d^2u}{dy^2} = 0")
    st.markdown("La ecuación se resuelve aplicando un esquema de diferencias finitas y/o aplicando el método del tiro (*shooting method*).")
    st.markdown("Los resultados numéricos se comparan con la solución analítica de la ecuación diferencial.")

    # Configuración de la barra lateral para la entrada de parámetros
    st.sidebar.header("Parámetros de cálculo")
    st.sidebar.subheader("Parámetros del flujo")
    uinf = st.sidebar.number_input("Velocidad de la corriente libre ($U_\infty$) en m/s", value=1.0, min_value=0.00,
                                   max_value=100.00, format="%.3f")
    vw = st.sidebar.number_input("Velocidad de succión ($V_w$) en m/s", value=0.1, min_value=0.001, max_value=10.000, 
                                 format="%.3f")
    mu = st.sidebar.number_input("Viscosidad absoluta ($\mu$) en $10^{-5}$ Pa·s", value=1.789, min_value=0.001, max_value=10.000, 
                                 format="%.3f") * 10**-5
    rho = st.sidebar.number_input("Densidad ($\rho$) en $\text{kg/m}^3$", value=1.225, min_value=0.001, max_value=10.000, 
                                  format="%.3f")
    nu = mu / rho

    st.sidebar.subheader("Parámetros de cálculo por diferencias finitas")
    n = st.sidebar.slider("Número de nodos", value=100, min_value=10, max_value=1000, step=10)
    esquema = st.sidebar.selectbox("Esquema de diferencias finitas", ["Hacia Atras de 1er orden", 
                "Hacia Adelante de 1er orden", "Central de 2do orden", "Hacia Atras de 2do orden", 
                "Hacia Adelante de 2do orden"])

    st.sidebar.subheader("Parámetros de cálculo por Shooting Method")
    a_inicial = st.sidebar.slider("Valor inicial de a", value=uinf * abs(vw) / nu, min_value=0.0, max_value=10.0, step=0.1)

    # Botón para ejecutar los cálculos
    if st.sidebar.button("🚀 Calcular Perfil de Velocidades"):
        with st.spinner("Resolviendo modelo numérico..."):
            # Definir el dominio espacial y el paso h
            delta_99 = 5.0 * (nu / vw)
            y_max = 2.0 * delta_99
            y = linspace(0, y_max, n)
            h = y[1] - y[0]

            # Ejecución de los métodos numéricos y analíticos
            udf = calcular_perfil_velocidades_df(uinf, vw, nu, n, h, esquema)
            usm = calcular_perfil_velocidades_sm(uinf, vw, nu, y, a_inicial)
            uan = calcular_perfil_velocidades_analitico(uinf, vw, nu, y)

            # Se guarda en el estado de la sesión para persistencia
            st.session_state.y = y
            st.session_state.udf = udf
            st.session_state.usm = usm
            st.session_state.uan = uan
            st.session_state.calculado = True

        st.success("¡Cálculo finalizado con éxito!")

    # Visualización de resultados (fuera del botón para mantenerlos en pantalla)
    if "calculado" in st.session_state and st.session_state.calculado:
        st.subheader("Resultados de la Simulación")

        fig, ax = subplots(figsize=(8, 5))
        ax.plot(
            st.session_state.udf,
            st.session_state.y,
            "b-o",
            label=f"Diferencias Finitas ({esquema})",
        )
        ax.plot(
            st.session_state.usm,
            st.session_state.y,
            "g-s",
            label="Método del Tiro (Shooting)",
        )
        ax.plot(
            st.session_state.uan,
            st.session_state.y,
            "r--",
            label="Solución Analítica",
        )
        ax.set_xlabel("Velocidad $u$ [m/s]")
        ax.set_ylabel("Coordenada $y$ [m]")
        ax.grid(True)
        ax.legend()

        st.pyplot(fig)

if __name__ == "__main__":
    main()
