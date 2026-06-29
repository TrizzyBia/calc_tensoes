import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc
import streamlit as st
import io
import contextlib

#===============================================================================================
#                                         ETAPA 1 
#===============================================================================================
#===============================================================================================
#                                         FUNÇÕES
#===============================================================================================

def criar_funcao_distribuida(expr):
    allowed = {
        "x": None,
        "np": np,
        "sin": np.sin,
        "cos": np.cos,
        "exp": np.exp,
        "sqrt": np.sqrt,
        "ln": np.log
    }

    def q(x):
        local_dict = allowed.copy()
        local_dict["x"] = x
        resultado = eval(expr, {"__builtins__": {}}, local_dict)

        # GARANTIA DE ARRAY
        if np.isscalar(resultado):
            return resultado * np.ones_like(x)

        return resultado

    return q

def diagramas_viga(dic,L,apoios):
    #dic --> list[dict]  (carregamentos)
    #L --> Float   (largura da viga)
    #apoios --> list[dict]  (tipos e posições dos apoios)
    
    #exemplo de dicionário força concentrada em y {"type": "point", "x": 2.0, "value": -10.0}

    #exemplo de dicionário momento concentrado {"type": "moment", "x": 2.0, "value": -10.0}

    #exemplo de dicionário carga distribuida {"type": "distributed", "x": [2.0, 3.0], "q": lambda x: 5+3x}
        #no exemplo da carga distribuida foi utilizada uma carga trapezoidal 5+3x

    #exemplo de dicionário força concentrada em x {"type": "point_x", "x": 2.0, "value": -10.0}

    #exemplo de dicionário carga distribuida {"type": "distributed_x", "x": [2.0, 3.0], "q": lambda x: 5+3x}
        #no exemplo da carga distribuida foi utilizada uma carga trapezoidal 5+3x

    #exemplo de dicionário apoios:
    #apoios = [
    #{"x": 0.0, "Rx": True,  "Ry": True,  "Mz": False},    pino
    #{"x": 4.0, "Rx": False, "Ry": True,  "Mz": False},    rolete
    #{"x": 2.0, "Rx": True, "Ry": True, "Mz": True}]    engaste

    
    # Checando se há apoios na mesma posição:
    posicao_apoios = []
    for a in apoios:
        posicao_apoios.append(round(a["x"],3)) #adiciona a posição de cada apoio em uma lista
    
    posicao_apoios2=list(set(posicao_apoios)) #define uma nova lista da posição dos apoios retirando os numeros repetidos
    
    if len(posicao_apoios2) < len(posicao_apoios):
        return("Há mais de um apoio no mesmo ponto, rever modelagem")

    # Checando a estaticidade da viga modelada
    n_Rx = sum(a["Rx"] for a in apoios)
    n_Ry = sum(a["Ry"] for a in apoios)
    n_Mz = sum(a["Mz"] for a in apoios)

    
    if n_Rx + n_Ry + n_Mz == 3:
        print("A viga é isostática")

    if n_Rx + n_Ry + n_Mz < 3 or n_Rx == 2 or n_Ry == 3 or n_Rx == 0:
        return("A viga modelada é hipoestática, reveja a modelagem dela")

    if n_Rx + n_Ry + n_Mz > 3: 
        return("A viga modelada é hiperestática, esse programa ainda não resolve esse tipo de viga")

    # Cálculo das reações de apoio:
    # F é o valor da força concentrada; r é a posição de cada força
    F_point=[]
    r_point=[]

    F_moment=[]
    r_moment=[]

    F_point_x=[]
    
    # no caso das cargas distribuidas, serão trazidas a informação das suas resultantes F_uniform_res 
    # e o ponto de aplicação delas no centroide r_uniform_res
    F_distributed_res=[]
    r_distributed_res=[]

    # Esse for serve para preencher o F e o r e testar se os inputs estão certos, percorrendo cada carregamento
    for load in dic:
        if load["type"]=="point":
            F_point.append(load["value"])
            r_point.append(load["x"])
        
        if load["type"]=="moment":
            F_moment.append(load["value"])
            r_moment.append(load["x"])
        
        if load["type"]=="distributed":
            x1=load["x"][0]
            x2=load["x"][1]
            x=np.linspace(x1,x2,1000)

            q = load["q"]

            #calculo da força resultante:
            
            #trazer os valores de q(x) para uma lista quando é uma distribuição que segue uma função
            if callable(q):
                qx = q(x)

            #trazer os valores de q(x) para uma lista quando é uma distribuição é constante
            else:
                qx = q * np.ones_like(x)

            #calculo da força resultante pela integral da carga
            F_distributed_res.append(-np.trapezoid(qx, x))

            #calculo do centroide para o caso da força resultante ser =0 ou diferente de 0
            if -np.trapezoid(qx, x)==0:
                r_distributed_res.append(0)
            
            else:
                r_distributed_res.append(np.trapezoid(qx*x, x)/np.trapezoid(qx, x))

        #forças em x para reação horizontal:
        if load["type"]=="point_x":
             F_point_x.append(load["value"])
        
        if load["type"]=="distributed_x":
            x1=load["x"][0]
            x2=load["x"][1]
            x=np.linspace(x1,x2,1000)
            
            q = load["q"]

            if callable(q):
                qx=q(x)     #Caso em que é uma carga distribuida variável
            
            else:
                qx = q * np.ones_like(x)  #Caso em que é uma carga distribuida constante

            F_point_x.append(np.trapezoid(qx,x)) #Força resultante em x da carga distribuida
        
        x = load["x"]

    # Caso: força/momento concentrado
        if np.isscalar(x):
            if x < 0 or x > L:
                return "Uma ou mais forças estão fora da viga, rever a modelagem do problema"

    # Caso: carga distribuída
        else:
            x1, x2 = x
            if x1 < 0 or x2 > L or x1 > x2:
                return "Uma ou mais forças estão fora da viga, rever a modelagem do problema"
         
    F_point=np.array(F_point,dtype=float)
    r_point=np.array(r_point,dtype=float)
    
    F_distributed_res=np.array(F_distributed_res,dtype=float)
    r_distributed_res=np.array(r_distributed_res,dtype=float)
    
    F_point_x=np.array(F_point_x,dtype=float)

    #Calculando Rx:
    Rx=-np.sum(F_point_x)
    for apoio in apoios:
        if apoio["Rx"] == True:
            xRx = apoio["x"]

    
    if n_Ry == 2: # Caso isostático: 2 reações em y
        x_apoio=[]
        for apoio in apoios:
            if apoio["Ry"] == True: 
                x_apoio.append(apoio["x"]) #Pegando o ponto de cada apoio em y
        x_apoio.sort() #Para organizar a ordem dos apoios
        
        #Calculando a reação em y no apoio da direita calculando momento em A:
        By=-(np.sum(F_point*(r_point-x_apoio[0]))+
             np.sum(F_moment)+
             np.sum(F_distributed_res*(r_distributed_res-x_apoio[0])))/(x_apoio[1]-x_apoio[0])

        #Calculando a reação em y de apoio na esquerda calculando somatório de forças em y:
        Ay=-(np.sum(F_point)+np.sum(F_distributed_res)) - By

        print(f"A reação em y do primeiro apoio = {Ay:.2f} kN")
                
        print(f"A reação em y do segundo apoio = {By:.2f} kN")
                
        print(f"A reação em x do apoio de 2º gênero = {Rx:.2f} kN")

    
    if len(apoios)==1: #engaste
        #Calculando a reação de apoio de momento no engaste calculando o momento em A(no apoio):
        x_apoio = apoios[0]["x"]
        Mz=-(np.sum(F_point*(r_point-x_apoio))+
             np.sum(F_moment)+
             np.sum(F_distributed_res*(r_distributed_res-x_apoio)))

        #Calculando a reação em y do engaste somatório de forças em y:
        Ay=-(np.sum(F_point)+np.sum(F_distributed_res))
        
        print(f"A reação em y = {Ay:.2f} kN")
        
        print(f"A reação em x = {Rx:.2f} kN")
        
        print(f"A reação de momento = {Mz:.2f} kNm")

        
    #Calculo dos esforços:
    #Definindo um eixo global:
    n = int(1000 * L)
    eixo=np.linspace(0,L,n)
    eixo[-1] = L
    
    #Definindo uma lista com a carga em todos o domínio da viga
    q_total_y = np.zeros_like(eixo)
    q_total_x = np.zeros_like(eixo)
    
    #Contabilizando a contribuição das cargas distribuidas para a carga total na viga
    for load in dic:
        if load["type"] == "distributed":
            x1, x2 = load["x"]
            dx = (eixo >= x1) & (eixo <= x2)
            q = load["q"]
            if callable(q):
                q_total_y[dx] += q(eixo[dx])
            else:
                q_total_y[dx] += q
        
        if load["type"] == "distributed_x":
            x1, x2 = load["x"]
            dx = (eixo >= x1) & (eixo <= x2)
            q = load["q"]
            if callable(q):
                q_total_x[dx] += q(eixo[dx])
            else:
                q_total_x[dx] += q
    
    
    #Calculando o cortante através da integral da carga em y
    Q = -sc.cumulative_trapezoid(q_total_y, eixo, initial=0)

    #Calculando o normal através da integral da carga em x
    N = -sc.cumulative_trapezoid(q_total_x, eixo, initial=0)
    
    #Adicionando a contribuição das reações de apoio no cortante
    if len(apoios) == 1:
        Q[eixo >= x_apoio] += Ay

    if n_Ry == 2:
        Q[eixo >= x_apoio[0]] += Ay
        Q[eixo >= x_apoio[1]] += By

    #Adicionando a contribuição das reações de apoio no normal
    N[eixo >= xRx] += -Rx
  
    #Adicionando a contribuição das forças concentradas no cortante e no normal
    for load in dic:
        if load["type"] == "point":
            Q[eixo >= load["x"]] += load["value"]
        
        if load["type"] == "point_x":
            N[eixo >= load["x"]] += -load["value"]
    
    #Calculando o momento fletor a partir da integral do cortante
    M = sc.cumulative_trapezoid(Q, eixo, initial=0)

    #Adicionando a contribuição dos momentos concentrados das cargas ou do engaste no momento fletor:
    for load in dic:
        if load["type"] == "moment":
            M[eixo >= load["x"]] += -load["value"]

    if len(apoios) == 1:
        M[eixo >= x_apoio] += -Mz
        

    #Preparando os dados de normal, cortanate e momento para que eles partam do 0 mesmo quando há descontinuidade gerada por reação ou carga concentrada:
    eixo_plot=np.insert(eixo,0,0)
    N_plot=np.insert(N,0,0)
    Q_plot=np.insert(Q,0,0)
    M_plot=np.insert(M,0,0)

    #Construção dos gráficos:
    fig, axs = plt.subplots(3, 1, figsize=(6, 6), sharex=True)
    fig.subplots_adjust(hspace=0.8)

    axs[0].plot(eixo_plot, N_plot)
    axs[0].set_title("Diagrama dos Esforços Normais")
    axs[0].set_ylabel("N (kN)")
    axs[0].axhline(0, color="black", linewidth=1)
    axs[0].grid()

    axs[1].plot(eixo_plot, Q_plot)
    axs[1].set_title("Diagrama dos Esforços Cortantes")
    axs[1].set_ylabel("Q (kN)")
    axs[1].axhline(0, color="black", linewidth=1)
    axs[1].grid()

    axs[2].plot(eixo_plot, -M_plot)
    axs[2].set_title("Diagrama dos Momentos Fletores")
    axs[2].set_xlabel("x (m)")
    axs[2].set_ylabel("M (kNm)")
    axs[2].axhline(0, color="black", linewidth=1)
    axs[2].grid()
    return { "fig" : fig, 
            "eixo" : eixo, 
            "N" : N, 
            "Q" : Q, 
            "M" : M, }

#===============================================================================================
#                                         INTERFACE
#===============================================================================================

if "etapa" not in st.session_state:
    st.session_state.etapa = 1

#Preparando as listas em que serão adicionados os apoios e os carregamentos
if "dic" not in st.session_state:
        st.session_state.dic = []

if "apoios" not in st.session_state:
        st.session_state.apoios = []

if st.session_state.etapa >= 1:

    st.title("Diagramas de Esforços em Vigas")

    st.write("Modelo simples de viga isostática")

    st.image("assets/eixos_adotados.jpeg", width=300)

    # Inputs
    L = st.number_input("Comprimento da viga (m)", min_value=0.1, value=5.0)


    st.subheader("Adicionar apoio")

    x_apoio = st.number_input("Posição do apoio (m)", min_value=0.0, value=0.0)

    tipo_apoio = st.selectbox(
        "Tipo de apoio",
        ["Rolete [Fy]", "Guia_horizontal[Fx]", "Pino[Fx,Fy]", "Engaste[Fx,Fy,Mz]"])

    if st.button("Adicionar apoio"):
        if x_apoio < 0 or x_apoio > L:
            st.error("O apoio está fora do domínio da viga!")
            st.stop()
        if tipo_apoio == "Pino[Fx,Fy]":
            apoio = {"x": x_apoio, "Rx": True, "Ry": True, "Mz": False}
        elif tipo_apoio == "Guia_horizontal[Fx]":
            apoio = {"x": x_apoio, "Rx": True, "Ry": False, "Mz": False}
        elif tipo_apoio == "Rolete [Fy]":
            apoio = {"x": x_apoio, "Rx": False, "Ry": True, "Mz": False}
        elif tipo_apoio == "Engaste[Fx,Fy,Mz]":
            apoio = {"x": x_apoio, "Rx": True, "Ry": True, "Mz": True}

        st.session_state.apoios.append(apoio)

    st.subheader("Apoios definidos")

    for i, apoio in enumerate(st.session_state.apoios):
        st.write(
            f"Apoio {i+1} → x = {apoio['x']} m | "
            f"Rx = {apoio['Rx']} | Ry = {apoio['Ry']} | Mz = {apoio['Mz']}"
        )

    if st.button("Limpar apoios"):
        st.session_state.apoios = []


    st.subheader("Carregamentos na viga")
    tipo_load = st.selectbox(
        "Tipo do carregamento",
        ["X pontual", "Y pontual", "Momento concentrado", "X distribuido", "Y distribuido"])

    if tipo_load in ["X distribuido", "Y distribuido"]:

        col1, col2 = st.columns(2)

        with col1:
            x_ini = st.number_input(
                "x inicial (m)",
                min_value=0.0,
                key="x_ini"
            )

        with col2:
            x_fim = st.number_input(
                "x final (m)",
                min_value=x_ini,
                key="x_fim"
            )

        x_load = [x_ini, x_fim]
        
        expr = st.text_input(
            "Função da carga q(x)",
            value="30*x",
            help="Exemplos: 10 | 30*x | 5 + 2*x | 20*sin(x)")   

    elif tipo_load in ["X pontual", "Y pontual"]:
        x_load = st.number_input(
            "Onde será aplicado o carregamento (m)",
            min_value=0.0,
            key="x_pontual")
        valor = st.number_input("Valor da força concentrada(kN))")  

    elif tipo_load == "Momento concentrado":
        x_load = st.number_input(
            "Onde será aplicado o carregamento (m)",
            min_value=0.0,
            key="x_pontual")
        valor = st.number_input("Valor do momento concentrado(kNm)")

    if st.button("Adicionar carregamento"):
        
        if tipo_load == "Y pontual":
            load = {"type": "point",
                    "x": x_load,
                    "value": valor}

        elif tipo_load == "Momento concentrado":
            load = {"type": "moment",
                    "x": x_load,
                    "value": valor}

        elif tipo_load == "X pontual":
            load = {"type": "point_x",
                    "x": x_load,
                    "value": valor}
            
        elif tipo_load == "Y distribuido":
            q_func = criar_funcao_distribuida(expr)

            # teste rápido de validade
            q_func(np.array([x_ini, x_fim]))

            load = {"type": "distributed",
                    "x": [x_ini, x_fim],
                    "q": q_func}

        elif tipo_load == "X distribuido":
            q_func = criar_funcao_distribuida(expr)

            q_func(np.array([x_ini, x_fim]))

            load = {"type": "distributed_x",
                    "x": [x_ini, x_fim],
                    "q": q_func}
        st.session_state.dic.append(load)

    st.subheader("Carregamentos adicionados")

    for i, load in enumerate(st.session_state.dic):
        st.write(f"Carregamento {i+1} → x = {load}")

    if st.button("Limpar carregamentos"):
        st.session_state.dic = []

    #Adptando returns para o streamlit:
    def rodar_viga_streamlit(dic, L, apoios):
        buffer = io.StringIO()

        try:
            with contextlib.redirect_stdout(buffer):
                resultado = diagramas_viga(dic, L, apoios)

            logs = buffer.getvalue()

            if isinstance(resultado, str):
                return {"status": "erro", "mensagem": resultado, "logs": logs}

            return {"status": "ok", 
                    "fig": resultado["fig"], 
                    "eixo": resultado["eixo"], 
                    "N": resultado["N"],
                    "Q": resultado["Q"], 
                    "M": resultado["M"], 
                    "logs": logs}

        except Exception as e:
            return {"status": "erro", "mensagem": str(e), "logs": buffer.getvalue()}
    
    #Preparando o resultado:
    if st.button("Calcular"):
        st.session_state.saida = rodar_viga_streamlit(
        st.session_state.dic,
        L,
        st.session_state.apoios)

    st.divider()

    if "saida" in st.session_state:
        saida = st.session_state.saida

        if saida["logs"]:
            st.subheader("📋 Resultados")
            st.text(saida["logs"])

        if saida["status"] == "erro":
            st.error(saida["mensagem"])
            st.stop()

        st.pyplot(saida["fig"])

        if st.button("Prosseguir para propriedades geométricas"):
            st.session_state.etapa = 2
            st.rerun()

#===============================================================================================
#                                         ETAPA 2 
#===============================================================================================

#===============================================================================================
#                                         FUNÇÕES
#===============================================================================================

def area_aba_inf(tfi, bfi):
    return tfi * bfi

def area_aba_sup(tfs, bfs):
    return tfs * bfs

def area_alma(h, tw):
    return h * tw

def area_total(h, tw, bfi, tfi, bfs, tfs):
    total = area_alma(h, tw) + area_aba_inf(tfi, bfi) + area_aba_sup(tfs, bfs)
    return total

# CENTROIDE
# Adotando o eixo y no eixo de simetria do perfil, temos que calcular yc
# Adotando o eixo z coincidindo com a parte inferior do perfil, temos que zc=0

# Calculando o yc:
def y_c(h, tw, bfi, tfi, bfs, tfs):
    yc=(area_aba_inf(tfi, bfi)*tfi/2+area_alma(h, tw)*(tfi+h/2)+area_aba_sup(tfs, bfs)*(tfi+h+tfs/2))/area_total(h, tw, bfi, tfi, bfs, tfs)
    return(yc)

# MOMENTO DE INÉRCIA
# eixo y vertical no centro de massa da viga e eixo z horizontal no centro de massa da viga

def mi_eixoy(h, tw, bfi, tfi, bfs, tfs):
    Iabay_inf = tfi * bfi**3 / 12
    Ialmay = h * tw**3 / 12
    Iabay_sup = tfs * bfs**3 / 12
    Itotaly = Iabay_inf+Ialmay+Iabay_sup
    return Itotaly

def mi_eixoz(h, tw, bfi, tfi, bfs, tfs, yc):
    d_abaz_inf = yc - tfi/2
    d_almaz = yc - (tfi + h/2)
    d_abaz_sup = yc - (tfi + h + tfs/2)
    
    Iabaz_inf = bfi * tfi**3 / 12 + area_aba_inf(tfi, bfi) * d_abaz_inf**2
    Ialmaz = tw * h**3 / 12 + area_alma(h, tw) * d_almaz**2
    Iabaz_sup = bfs * tfs**3 / 12 + area_aba_sup(tfs, bfs) * d_abaz_sup**2

    Itotalz = Iabaz_inf + Ialmaz + Iabaz_sup
   
    return Itotalz

#===============================================================================================
#                                         INTERFACE  
#===============================================================================================

if st.session_state.etapa >= 2:
    st.header("Propriedades Geométricas")

    col1, col2 = st.columns([3,2])

    with col1:
        bfs = st.number_input('Largura da aba superior - bfs (mm)', min_value=1.0)
        tfs = st.number_input('Espessura da aba superior - tfs (mm)', min_value=1.0)
        h = st.number_input('Altura interna da alma - h (mm)', min_value=1.0)
        tw = st.number_input('Espessura da alma - tw (mm)', min_value=1.0)
        bfi = st.number_input('Largura da aba inferior - bfi (mm)', min_value=1.0)
        tfi = st.number_input('Espessura da aba inferior - tfi (mm)', min_value=1.0)

    with col2:
        st.markdown("<h4 style='text-align: center;'>Geometria da Seção</h4>", unsafe_allow_html=True)

        st.image("assets/perfil_i.png.png", use_container_width=True)

    if st.button('Calcular Propriedades'):
        # salvando para a próxima etapa
        st.session_state.area = area_total(h, tw, bfi, tfi, bfs, tfs)
        st.session_state.yc = y_c(h, tw, bfi, tfi, bfs, tfs)
        st.session_state.miz = mi_eixoz(h, tw, bfi, tfi, bfs, tfs, st.session_state.yc)
        st.session_state.miy = mi_eixoy(h, tw, bfi, tfi, bfs, tfs)
        st.session_state.h = h
        st.session_state.tfi = tfi
        st.session_state.tw = tw
        st.session_state.bfi = bfi
        st.session_state.bfs = bfs
        st.session_state.tfs = tfs

    st.divider()

    if "area" in st.session_state:
        st.subheader('📋 Resultados')

        with st.container(border=True):
            st.subheader("Área da Seção Transversal")
            st.write(f"**A =** {st.session_state.area:.3f} mm²")

        with st.container(border=True):
            st.subheader("Coordenada y do Centróide")
            st.caption("***Em relação à base do perfil**")
            st.write(f"**yc =** {st.session_state.yc:.3f} mm")

        with st.container(border=True):
            st.subheader("Momento de Inércia (eixo z)")
            st.write(f"**Iz =** {st.session_state.miz:.3f} mm⁴")

        with st.container(border=True):
            st.subheader("Momento de Inércia (eixo y)")
            st.write(f"**Iy =** {st.session_state.miy:.3f} mm⁴")

        if st.button("Prosseguir para cálculo de tensões"):
            st.session_state.etapa = 3
            st.rerun()

#===============================================================================================
#                                         ETAPA 3 
#===============================================================================================
#===============================================================================================
#                                         FUNÇÕES
#===============================================================================================

# Tensões normais
def sigma_x(N, M, Iz, y, yc, area):
    normal = (N * 1e3) / area
    fletor = (M * 1e6 * (yc - y)) / Iz
    return normal + fletor

# Tensões cisalhantes

# Cálculo do Q
def q_estatico(y, yc, tfi, bfi, h, tw, tfs, bfs):
    if y <= yc:
        area_q = (tfi * bfi) + (tw * (y - tfi))
        yc_q = ((tfi * bfi * tfi / 2) + (tw * (y - tfi) * (tfi + (y - tfi) / 2))) / area_q
        return area_q * (yc - yc_q) 
    
    else:
        area_q = (tfs * bfs) + tw * (h + tfi - y)
        yc_q = (tfs * bfs*((tfi + h + tfs/2)) + (tw * (h + tfi - y) * (y + (h - y + tfi)/2)))/area_q
        return area_q * (yc_q - yc) 
    
def tau_xy(Q, Qest, Iz, tw):
    cortante = (Q * 1e3 * Qest)/(Iz * tw)
    return cortante

def tensoes_principais(sigma, tau):
    taumax = np.sqrt((sigma / 2)**2 + (tau**2))
    sigma1 = sigma/2 + taumax
    sigma2 = sigma/2 - taumax
    return sigma1, sigma2, taumax

#===============================================================================================
#                                         INTERFACE
#===============================================================================================

if st.session_state.etapa >= 3:
    # trazendo da etapa anterior
    area = st.session_state.area
    Iz = st.session_state.miz
    Iy = st.session_state.miy
    yc = st.session_state.yc
    h = st.session_state.h
    tw = st.session_state.tw
    tfi = st.session_state.tfi
    bfi = st.session_state.bfi
    tfs = st.session_state.tfs
    bfs = st.session_state.bfs

    st.header("Cálculo das Tensões")

    # ponto x
    x_consulta = st.number_input("Ponto x da viga (m)", min_value=0.0, max_value=L)
    
    if "saida" in st.session_state:
        eixo = st.session_state.saida["eixo"]
        N = st.session_state.saida["N"]
        Q = st.session_state.saida["Q"]
        M = st.session_state.saida["M"]

        descontinuidades = []

        # cargas
        for carga in st.session_state.dic:
            if carga["type"] == "point":
                descontinuidades.append(carga["x"])

            if carga["type"] == "moment":
                descontinuidades.append(carga["x"])

        # apoios
        for apoio in st.session_state.apoios:
            if apoio["Ry"] != 0 or apoio["Mz"] != 0:
                descontinuidades.append(apoio["x"])

        descontinuidades = sorted(set(descontinuidades))

        tol = 0.001  # 1 mm

        # verifica se o ponto é válido
        ponto_valido = not any(abs(x_consulta - x) < tol for x in descontinuidades)

        if not ponto_valido:
            # apaga resultados anteriores
            for chave in ("sigma", "tau", "sigma1", "sigma2", "taumax",):
                st.session_state.pop(chave, None)

            # avisa ao usuário sobre pontos inválidos
            st.warning("Este é um ponto de descontinuidade, escolha outro ponto.\n\n" 
                       f"Descontinuidades em: {', '.join(f'{x:.3f}' for x in descontinuidades)} m.")


        Nx = np.interp(x_consulta, eixo, N)
        Qx = np.interp(x_consulta, eixo, Q)
        Mx = np.interp(x_consulta, eixo, M)

        # ponto y 
        if ponto_valido:
            y = st.number_input("Ponto y da seção (mm)", min_value=tfi, max_value=tfi + h, value=tfi)
            st.caption("***Em relação à base do perfil**")

            if st.button("Calcular Tensões"):
                Qest = q_estatico(y, yc, tfi, bfi, h, tw, tfs, bfs)
                st.session_state.sigma = sigma_x(Nx, Mx, Iz, y, yc, area) 
                st.session_state.tau = tau_xy(Qx, Qest, Iz, tw)

                sigma1, sigma2, taumax = tensoes_principais(st.session_state.sigma, st.session_state.tau)
                st.session_state.sigma1 = sigma1
                st.session_state.sigma2 = sigma2
                st.session_state.taumax = taumax

    st.divider()

    if "sigma" in st.session_state:
        st.subheader("📋 Resultados")

        with st.container(border=True):
            st.subheader("Esforços na Seção")
            st.write(f"**N =** {Nx:.3f} kN")
            st.write(f"**V =** {Qx:.3f} kN")
            st.write(f"**M =** {Mx:.3f} kNm")

        with st.container(border=True):
            st.subheader("Tensões no Ponto Escolhido")
            st.write(f"**σx =** {st.session_state.sigma:.3f} MPa")
            st.write(f"**τxy =** {st.session_state.tau:.3f} MPa")


        with st.container(border=True):
            st.subheader("Tensões Principais")
            st.write(f"**σ₁ =** {st.session_state.sigma1:.3f} MPa")
            st.write(f"**σ₂ =** {st.session_state.sigma2:.3f} MPa")

        with st.container(border=True):
            st.subheader("Tensão de Cisalhamento Máxima")
            st.write(f"**τmáx =** {st.session_state.taumax:.3f} MPa")

        st.success("Tensões calculadas com sucesso!")
