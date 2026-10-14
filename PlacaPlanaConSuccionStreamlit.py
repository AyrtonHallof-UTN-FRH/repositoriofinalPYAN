import streamlit as st
from numpy import linalg, array, zeros
import matplotlib.pyplot as plt

#Configuracion de la pagina#
st.set_page_config(page_title="Perfil de Velocidades en una Placa Plana con Succión", 
                   page_icon="https://share.google/K4OvrbntlsuX37Y8H.png",
                   layout="wide",
                   initial_sidebar_state="expanded")

#Presentacion de la pagina#
st.title("Perfil de Velocidades en una Placa Plana con Succión")
st.subheader("Cátedra de Programación y Análisis Numérico-UTN FRH")
st.markdown("En esta aplicacion se puede obtener el perfil de velociades sobre una placa plana con succión,")
st.markdown("que viene determinado por la siguiente ecuación diferencial:")         
st.markdown("$$V_w \frac{du}{dy} + \nu \frac{d^2u}{dy^2} = 0$$")
st.markdown("la ecuacion se resuleve aplicando un esquema de diferencias finitas y/o aplicando (*shooting method*).")
st.markdown("Los resultados numericos se comparan con la solucion analitica de la ecuacion diferencial")

#configuracion de la barra lateral#
st.sidebar.header("Parametros de cálculo")
st.sidebar.subheader("Parametros del flujo")
uinf= st.sidebar.number_input("Velocidad de la corriente libre $$U_\infty$$ en m/s", value=1.0, min_value=0.00, max_value=200.00, format="%.3f")
vw= st.sidebar.number_input("Velocidad de succión $$V_w$$ en m/s", value=0.1, min_value=0.001, max_value=10.000, format="%.3f")
mu= st.sidebar.number_input("Viscosidad absoluta $$\mu$$ en 10^-5 Pa·s", value=1.789, min_value=0.001, max_value=10.000, format="%.3f")*10**-5
rho= st.sidebar.number_input("Densidad $$\rho$$ en kg/m^3", value=1.225, min_value=0.001, max_value=10.000, format="%.3f")
nu= mu/rho

st.sidebar.subheader("Parametros de cálculo numérico")
n= st.sidebar.slider("Número de nodos", value=100, min_value=10, max_value=1000, step=10)
esquema= st.sidebar.selectbox("Esquema de diferencias finitas", ["Hacia Atras de 1er orden", "Hacia Adelante de 1er orden", "Central de 2do orden", "Hacia Atras de 2do orden", "Hacia Adelante de 2do orden"])