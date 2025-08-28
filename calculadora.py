import tkinter as tk
from tkinter import ttk, messagebox
import math
import statistics
try:
    import numpy as np
    NUMPY_DISPONIVEL = True
except ImportError:
    NUMPY_DISPONIVEL = False
    messagebox.showwarning("Aviso", "Biblioteca NumPy não encontrada. Algumas funcionalidades estarão limitadas.")

# --------- Funções Matemáticas Gerais ---------
def bhaskara(a, b, c):
    delta = b*b - 4*a*c
    if delta < 0:
        return delta, None, None
    x1 = (-b + math.sqrt(delta)) / (2*a)
    x2 = (-b - math.sqrt(delta)) / (2*a)
    return delta, x1, x2

def area_circulo(r):
    return math.pi * r*r

def circunferencia(r):
    return 2 * math.pi * r

def diametro(r):
    return 2 * r

def pythagoras_hip(cat1, cat2):
    return math.sqrt(cat1*cat1 + cat2*cat2)

def pythagoras_cat(hip, cat1):
    if hip <= cat1:
        return None
    return math.sqrt(hip*hip - cat1*cat1)

def juros_compostos(capital, taxa, n):
    montante = capital * (1 + taxa)**n
    juros = montante - capital
    return juros, montante

# --------- Funções de Geometria ---------
def area_poligono_regular(n_lados, comprimento_lado):
    return (n_lados * comprimento_lado**2) / (4 * math.tan(math.pi / n_lados))

def angulo_interno_poligono(n_lados):
    return (n_lados - 2) * 180 / n_lados

def angulo_externo_poligono(n_lados):
    return 360 / n_lados

def area_superficie_prisma_retangular(l, w, h):
    return 2 * (l*w + l*h + w*h)

def volume_prisma_retangular(l, w, h):
    return l * w * h

def area_superficie_prisma_triangular(a, b, c, h):
    s = (a + b + c) / 2
    area_base = math.sqrt(s * (s - a) * (s - b) * (s - c))
    area_lateral = (a + b + c) * h
    return area_lateral + 2 * area_base

def volume_prisma_triangular(a, b, c, h):
    s = (a + b + c) / 2
    area_base = math.sqrt(s * (s - a) * (s - b) * (s - c))
    return area_base * h

def area_superficie_cilindro(r, h):
    return 2 * math.pi * r * (r + h)

def volume_cilindro(r, h):
    return math.pi * r * r * h

def area_superficie_esfera(r):
    return 4 * math.pi * r * r

def volume_esfera(r):
    return (4/3) * math.pi * r * r * r

def area_superficie_cone(r, h):
    return math.pi * r * (r + math.sqrt(h*h + r*r))

def volume_cone(r, h):
    return (1/3) * math.pi * r * r * h

# --------- Funções de Trigonometria ---------
def lei_dos_senos(a, A, b, B, c, C):
    # Resolve triângulo usando lei dos senos
    if a and A:
        if b and not B:
            B = math.degrees(math.asin((b * math.sin(math.radians(A))) / a))
        if c and not C:
            C = math.degrees(math.asin((c * math.sin(math.radians(A))) / a))
    return A, B, C, a, b, c

def lei_dos_cossenos(a, b, c, A, B, C):
    # Resolve triângulo usando lei dos cossenos
    if a and b and c and not A:
        A = math.degrees(math.acos((b*b + c*c - a*a) / (2*b*c)))
    if a and b and C and not c:
        c = math.sqrt(a*a + b*b - 2*a*b*math.cos(math.radians(C)))
    return A, B, C, a, b, c

# --------- Funções de Álgebra ---------
def progressao_aritmetica(a1, n, r):
    an = a1 + (n - 1) * r
    sn = n * (a1 + an) / 2
    return an, sn

def progressao_geometrica(a1, n, q):
    an = a1 * (q ** (n - 1))
    if q == 1:
        sn = n * a1
    else:
        sn = a1 * (1 - q**n) / (1 - q)
    return an, sn

def fatorial(n):
    if n == 0:
        return 1
    return n * fatorial(n-1)

def combinacao(n, k):
    return fatorial(n) / (fatorial(k) * fatorial(n - k))

def permutacao(n, k=None):
    if k is None:
        return fatorial(n)
    return fatorial(n) / fatorial(n - k)

def arranjo(n, k):
    return fatorial(n) / fatorial(n - k)

# --------- Funções de Probabilidade ---------
def probabilidade_evento(casos_favoraveis, casos_possiveis):
    return casos_favoraveis / casos_possiveis

def probabilidade_condicional(p_a_e_b, p_b):
    return p_a_e_b / p_b

def probabilidade_bayes(p_b_dado_a, p_a, p_b):
    return (p_b_dado_a * p_a) / p_b

# --------- Funções de Matrizes (se numpy disponível) ---------
if NUMPY_DISPONIVEL:
    def matriz_soma(matriz_a, matriz_b):
        return np.add(matriz_a, matriz_b)

    def matriz_subtracao(matriz_a, matriz_b):
        return np.subtract(matriz_a, matriz_b)

    def matriz_multiplicacao(matriz_a, matriz_b):
        return np.matmul(matriz_a, matriz_b)

    def matriz_transposta(matriz):
        return np.transpose(matriz)

    def matriz_determinante(matriz):
        return np.linalg.det(matriz)

    def matriz_inversa(matriz):
        try:
            return np.linalg.inv(matriz)
        except:
            return None  # Matriz singular

    # --------- Funções de Sistemas Lineares ---------
    def resolver_sistema_lineares(matriz_coeficientes, vetor_constantes):
        try:
            return np.linalg.solve(matriz_coeficientes, vetor_constantes)
        except:
            return None  # Sistema sem solução única
else:
    # Funções placeholder quando numpy não está disponível
    def matriz_soma(matriz_a, matriz_b):
        return "Numpy não disponível"
    
    def matriz_subtracao(matriz_a, matriz_b):
        return "Numpy não disponível"
    
    def matriz_multiplicacao(matriz_a, matriz_b):
        return "Numpy não disponível"
    
    def matriz_transposta(matriz):
        return "Numpy não disponível"
    
    def matriz_determinante(matriz):
        return "Numpy não disponível"
    
    def matriz_inversa(matriz):
        return "Numpy não disponível"
    
    def resolver_sistema_lineares(matriz_coeficientes, vetor_constantes):
        return "Numpy não disponível"

# --------- Funções de Cálculo ---------
def derivada_polinomial(coeficientes):
    # coeficientes: lista [a, b, c, ...] para ax^n + bx^(n-1) + ...
    derivada = []
    grau = len(coeficientes) - 1
    for i, coef in enumerate(coeficientes):
        if i < grau:
            derivada.append(coef * (grau - i))
    return derivada

def integral_polinomial(coeficientes):
    # coeficientes: lista [a, b, c, ...] para ax^n + bx^(n-1) + ...
    integral = []
    grau = len(coeficientes) - 1
    for i, coef in enumerate(coeficientes):
        integral.append(coef / (grau - i + 1))
    integral.append(0)  # constante de integração
    return integral

# --------- Funções de Estatística ---------
def medidas_estatisticas(dados):
    dados = [float(x) for x in dados.split()]
    media = statistics.mean(dados)
    mediana = statistics.median(dados)
    try:
        moda = statistics.mode(dados) if len(dados) > 1 and max([dados.count(x) for x in dados]) > 1 else "Nenhuma"
    except:
        moda = "Nenhuma"
    variancia = statistics.variance(dados)
    desvio_padrao = statistics.stdev(dados)
    return media, mediana, moda, variancia, desvio_padrao

# --------- Funções de Matemática Financeira ---------
def juros_simples(capital, taxa, tempo):
    juros = capital * taxa * tempo
    montante = capital + juros
    return juros, montante

def valor_prestacao_parcela(valor, taxa, periodos):
    if taxa == 0:
        return valor / periodos
    return valor * (taxa * (1 + taxa)**periodos) / ((1 + taxa)**periodos - 1)

# --------- Funções do Teclado ---------
def inserir_texto(texto):
    widget_focado = janela.focus_get()
    if isinstance(widget_focado, tk.Entry):
        posicao = widget_focado.index(tk.INSERT)
        widget_focado.insert(posicao, texto)
        atualizar_previa()

def apagar_caractere():
    widget_focado = janela.focus_get()
    if isinstance(widget_focado, tk.Entry):
        posicao = widget_focado.index(tk.INSERT)
        if posicao > 0:
            widget_focado.delete(posicao-1, posicao)
            atualizar_previa()

def limpar_campo():
    widget_focado = janela.focus_get()
    if isinstance(widget_focado, tk.Entry):
        widget_focado.delete(0, tk.END)
        atualizar_previa()

def limpar_tudo():
    for entry in [entry1, entry2, entry3, entry4, entry5, entry6, entry_pyth1, entry_pyth2]:
        entry.delete(0, tk.END)
    resultado.set("")
    preview.set("")
    atualizar_previa()

def inserir_operador(operador):
    widget_focado = janela.focus_get()
    if isinstance(widget_focado, tk.Entry):
        posicao = widget_focado.index(tk.INSERT)
        widget_focado.insert(posicao, operador)
        atualizar_previa()

def calcular_expressao(expressao):
    try:
        expressao = expressao.replace(",", ".").replace("^", "**")
        return eval(expressao)
    except:
        return None

def atualizar_previa(event=None):
    if combo.get() == "Calculadora normal":
        expressao = entry1.get()
        if expressao:
            res = calcular_expressao(expressao)
            if res is not None:
                preview.set(f"= {res}")
            else:
                preview.set("")
        else:
            preview.set("")

# --------- Atualiza rótulos conforme a operação ---------
def atualizar_labels(event=None):
    opcao = combo.get()
    
    # Esconder todos os widgets primeiro (com tratamento de erro)
    widgets_para_esconder = [label2, entry2, label3, entry3, label4, entry4, label5, entry5, label6, entry6,
                            op_circ, op_pyth, op_poligono, op_prisma, op_trig, op_alg, op_fin, op_prob, op_matriz]
    
    for widget in widgets_para_esconder:
        try:
            widget.pack_forget()
        except:
            pass  # Ignora widgets que não existem ou não estão empacotados
    
    # Reposicionar label1 e entry1
    label1.pack(pady=2)
    entry1.pack(pady=2)
    
    # Configurar conforme a operação selecionada
    if opcao == "Equação 2º Grau":
        label1.config(text="a:")
        label2.config(text="b:")
        label3.config(text="c:")
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)

    elif opcao == "Círculo":
        label1.config(text="Raio:")
        op_circ.pack(pady=5)

    elif opcao == "Pitágoras":
        label1.pack_forget()
        entry1.pack_forget()
        op_pyth.pack(pady=5)
        atualizar_pitagoras()

    elif opcao == "Juros compostos":
        label1.config(text="Capital (C):")
        label2.config(text="Taxa (i) decimal:")
        label3.config(text="Períodos (n):")
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)

    elif opcao == "Polígonos Regulares":
        op_poligono.pack(pady=5)
        atualizar_poligono()

    elif opcao == "Prismas/Sólidos":
        op_prisma.pack(pady=5)
        atualizar_prisma()

    elif opcao == "Trigonometria":
        op_trig.pack(pady=5)
        atualizar_trigonometria()

    elif opcao == "Álgebra":
        op_alg.pack(pady=5)
        atualizar_algebra()

    elif opcao == "Probabilidade":
        op_prob.pack(pady=5)
        atualizar_probabilidade()

    elif opcao == "Matrizes":
        op_matriz.pack(pady=5)
        atualizar_matrizes()

    elif opcao == "Sistemas Lineares":
        label1.config(text="Número de variáveis:")
        entry1.pack(pady=2)

    elif opcao == "Estatística":
        label1.config(text="Dados (separados por espaço):")
        entry1.pack(pady=2)

    elif opcao == "Matemática Financeira":
        op_fin.pack(pady=5)
        atualizar_financeira()

    elif opcao == "Cálculo":
        label1.config(text="Coeficientes (separados por espaço):")
        entry1.pack(pady=2)

    elif opcao == "Calculadora normal":
        label1.config(text="Digite a expressão:")
        entry1.delete(0, tk.END)
    
    # Limpar preview quando mudar de operação
    preview.set("")
    resultado.set("")

def atualizar_pitagoras():
    escolha = op_pyth.get()
    label_pyth1.pack_forget()
    entry_pyth1.pack_forget()
    label_pyth2.pack_forget()
    entry_pyth2.pack_forget()
    
    if escolha == "Hipotenusa":
        label_pyth1.config(text="Cateto 1:")
        label_pyth2.config(text="Cateto 2:")
    else:
        label_pyth1.config(text="Hipotenusa:")
        label_pyth2.config(text="Cateto conhecido:")
    
    label_pyth1.pack(pady=2)
    entry_pyth1.pack(pady=2)
    label_pyth2.pack(pady=2)
    entry_pyth2.pack(pady=2)

def atualizar_poligono():
    escolha = op_poligono.get()
    
    for widget in [label1, entry1, label2, entry2, label3, entry3]:
        widget.pack_forget()
    
    if escolha == "Área":
        label1.config(text="Número de lados:")
        label2.config(text="Comprimento do lado:")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
    elif escolha == "Ângulos":
        label1.config(text="Número de lados:")
        label1.pack(pady=2)
        entry1.pack(pady=2)

def atualizar_prisma():
    escolha = op_prisma.get()
    
    for widget in [label1, entry1, label2, entry2, label3, entry3, label4, entry4]:
        widget.pack_forget()
    
    if escolha == "Prisma Retangular":
        label1.config(text="Comprimento (l):")
        label2.config(text="Largura (w):")
        label3.config(text="Altura (h):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)
    elif escolha == "Prisma Triangular":
        label1.config(text="Lado A da base:")
        label2.config(text="Lado B da base:")
        label3.config(text="Lado C da base:")
        label4.config(text="Altura (h):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)
        label4.pack(pady=2)
        entry4.pack(pady=2)
    elif escolha == "Cilindro":
        label1.config(text="Raio (r):")
        label2.config(text="Altura (h):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
    elif escolha == "Esfera":
        label1.config(text="Raio (r):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
    elif escolha == "Cone":
        label1.config(text="Raio (r):")
        label2.config(text="Altura (h):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)

def atualizar_trigonometria():
    escolha = op_trig.get()
    
    for widget in [label1, entry1, label2, entry2, label3, entry3, label4, entry4, label5, entry5, label6, entry6]:
        widget.pack_forget()
    
    if escolha == "Lei dos Senos":
        label1.config(text="Lado a:")
        label2.config(text="Ângulo A:")
        label3.config(text="Lado b:")
        label4.config(text="Ângulo B:")
        label5.config(text="Lado c:")
        label6.config(text="Ângulo C:")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)
        label4.pack(pady=2)
        entry4.pack(pady=2)
        label5.pack(pady=2)
        entry5.pack(pady=2)
        label6.pack(pady=2)
        entry6.pack(pady=2)
    elif escolha == "Lei dos Cossenos":
        label1.config(text="Lado a:")
        label2.config(text="Lado b:")
        label3.config(text="Lado c:")
        label4.config(text="Ângulo A:")
        label5.config(text="Ângulo B:")
        label6.config(text="Ângulo C:")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)
        label4.pack(pady=2)
        entry4.pack(pady=2)
        label5.pack(pady=2)
        entry5.pack(pady=2)
        label6.pack(pady=2)
        entry6.pack(pady=2)

def atualizar_algebra():
    escolha = op_alg.get()
    
    for widget in [label1, entry1, label2, entry2, label3, entry3]:
        widget.pack_forget()
    
    if escolha == "Progressão Aritmética":
        label1.config(text="Primeiro termo (a1):")
        label2.config(text="Número de termos (n):")
        label3.config(text="Razão (r):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)
    elif escolha == "Progressão Geométrica":
        label1.config(text="Primeiro termo (a1):")
        label2.config(text="Número de termos (n):")
        label3.config(text="Razão (q):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)
    elif escolha == "Análise Combinatória":
        label1.config(text="n:")
        label2.config(text="k:")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)

def atualizar_probabilidade():
    escolha = op_prob.get()
    
    for widget in [label1, entry1, label2, entry2, label3, entry3]:
        widget.pack_forget()
    
    if escolha == "Probabilidade Simples":
        label1.config(text="Casos favoráveis:")
        label2.config(text="Casos possíveis:")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
    elif escolha == "Probabilidade Condicional":
        label1.config(text="P(A e B):")
        label2.config(text="P(B):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
    elif escolha == "Teorema de Bayes":
        label1.config(text="P(B|A):")
        label2.config(text="P(A):")
        label3.config(text="P(B):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)

def atualizar_matrizes():
    escolha = op_matriz.get()
    
    for widget in [label1, entry1, label2, entry2, label3, entry3]:
        widget.pack_forget()
    
    if escolha in ["Soma", "Subtração", "Multiplicação"]:
        label1.config(text="Matriz A (linhas separadas por ;):")
        label2.config(text="Matriz B (linhas separadas por ;):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
    elif escolha in ["Transposta", "Determinante", "Inversa"]:
        label1.config(text="Matriz (linhas separadas por ;):")
        label1.pack(pady=2)
        entry1.pack(pady=2)

def atualizar_financeira():
    escolha = op_fin.get()
    
    for widget in [label1, entry1, label2, entry2, label3, entry3]:
        widget.pack_forget()
    
    if escolha == "Juros Simples":
        label1.config(text="Capital (C):")
        label2.config(text="Taxa (i) decimal:")
        label3.config(text="Tempo (n):")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)
    elif escolha == "Prestação Fixa":
        label1.config(text="Valor financiado:")
        label2.config(text="Taxa (i) decimal:")
        label3.config(text="Número de parcelas:")
        label1.pack(pady=2)
        entry1.pack(pady=2)
        label2.pack(pady=2)
        entry2.pack(pady=2)
        label3.pack(pady=2)
        entry3.pack(pady=2)

# --------- Função para cálculo ---------
def calcular():
    try:
        opcao = combo.get()

        if opcao == "Equação 2º Grau":
            a = float(entry1.get().replace(",", "."))
            b = float(entry2.get().replace(",", "."))
            c = float(entry3.get().replace(",", "."))
            delta, x1, x2 = bhaskara(a, b, c)
            if x1 is None:
                resultado.set(f"Δ = {delta}\nSem raízes reais.")
            else:
                resultado.set(f"Δ = {delta}\nRaízes:\nx1 = {x1:.4f}, x2 = {x2:.4f}")

        elif opcao == "Círculo":
            r = float(entry1.get().replace(",", "."))
            escolha = op_circ.get()
            if escolha == "Área":
                res = area_circulo(r)
            elif escolha == "Circunferência":
                res = circunferencia(r)
            elif escolha == "Diámetro":
                res = diametro(r)
            resultado.set(f"{escolha} = {res:.4f}")

        elif opcao == "Pitágoras":
            escolha = op_pyth.get()
            if escolha == "Hipotenusa":
                c1 = float(entry_pyth1.get().replace(",", "."))
                c2 = float(entry_pyth2.get().replace(",", "."))
                res = pythagoras_hip(c1, c2)
                resultado.set(f"Hipotenusa = {res:.4f}")
            elif escolha == "Cateto faltante":
                hip = float(entry_pyth1.get().replace(",", "."))
                c1 = float(entry_pyth2.get().replace(",", "."))
                res = pythagoras_cat(hip, c1)
                if res is None:
                    resultado.set("Erro: Hipotenusa deve ser maior que o cateto.")
                else:
                    resultado.set(f"Cateto faltante = {res:.4f}")

        elif opcao == "Juros compostos":
            C = float(entry1.get().replace(",", "."))
            i = float(entry2.get().replace(",", "."))
            n = int(entry3.get())
            j, m = juros_compostos(C, i, n)
            resultado.set(f"Juros = R$ {j:.2f}\nMontante = R$ {m:.2f}")

        elif opcao == "Polígonos Regulares":
            escolha = op_poligono.get()
            if escolha == "Área":
                n = int(entry1.get())
                lado = float(entry2.get().replace(",", "."))
                area = area_poligono_regular(n, lado)
                ang_int = angulo_interno_poligono(n)
                ang_ext = angulo_externo_poligono(n)
                resultado.set(f"Área = {area:.4f}\nÂngulo interno = {ang_int:.2f}°\nÂngulo externo = {ang_ext:.2f}°")
            elif escolha == "Ângulos":
                n = int(entry1.get())
                ang_int = angulo_interno_poligono(n)
                ang_ext = angulo_externo_poligono(n)
                resultado.set(f"Ângulo interno = {ang_int:.2f}°\nÂngulo externo = {ang_ext:.2f}°")

        elif opcao == "Prismas/Sólidos":
            escolha = op_prisma.get()
            if escolha == "Prisma Retangular":
                l = float(entry1.get().replace(",", "."))
                w = float(entry2.get().replace(",", "."))
                h = float(entry3.get().replace(",", "."))
                area = area_superficie_prisma_retangular(l, w, h)
                volume = volume_prisma_retangular(l, w, h)
                resultado.set(f"Área da superfície = {area:.4f}\nVolume = {volume:.4f}")
            elif escolha == "Prisma Triangular":
                a = float(entry1.get().replace(",", "."))
                b = float(entry2.get().replace(",", "."))
                c = float(entry3.get().replace(",", "."))
                h = float(entry4.get().replace(",", "."))
                area = area_superficie_prisma_triangular(a, b, c, h)
                volume = volume_prisma_triangular(a, b, c, h)
                resultado.set(f"Área da superfície = {area:.4f}\nVolume = {volume:.4f}")
            elif escolha == "Cilindro":
                r = float(entry1.get().replace(",", "."))
                h = float(entry2.get().replace(",", "."))
                area = area_superficie_cilindro(r, h)
                volume = volume_cilindro(r, h)
                resultado.set(f"Área da superfície = {area:.4f}\nVolume = {volume:.4f}")
            elif escolha == "Esfera":
                r = float(entry1.get().replace(",", "."))
                area = area_superficie_esfera(r)
                volume = volume_esfera(r)
                resultado.set(f"Área da superfície = {area:.4f}\nVolume = {volume:.4f}")
            elif escolha == "Cone":
                r = float(entry1.get().replace(",", "."))
                h = float(entry2.get().replace(",", "."))
                area = area_superficie_cone(r, h)
                volume = volume_cone(r, h)
                resultado.set(f"Área da superfície = {area:.4f}\nVolume = {volume:.4f}")

        elif opcao == "Trigonometria":
            escolha = op_trig.get()
            if escolha == "Lei dos Senos":
                a = float(entry1.get().replace(",", ".")) if entry1.get() else None
                A = float(entry2.get().replace(",", ".")) if entry2.get() else None
                b = float(entry3.get().replace(",", ".")) if entry3.get() else None
                B = float(entry4.get().replace(",", ".")) if entry4.get() else None
                c = float(entry5.get().replace(",", ".")) if entry5.get() else None
                C = float(entry6.get().replace(",", ".")) if entry6.get() else None
                
                A, B, C, a, b, c = lei_dos_senos(a, A, b, B, c, C)
                resultado.set(f"Ângulo A = {A:.2f}°\nÂngulo B = {B:.2f}°\nÂngulo C = {C:.2f}°\nLado a = {a:.4f}\nLado b = {b:.4f}\nLado c = {c:.4f}")
                
            elif escolha == "Lei dos Cossenos":
                a = float(entry1.get().replace(",", ".")) if entry1.get() else None
                b = float(entry2.get().replace(",", ".")) if entry2.get() else None
                c = float(entry3.get().replace(",", ".")) if entry3.get() else None
                A = float(entry4.get().replace(",", ".")) if entry4.get() else None
                B = float(entry5.get().replace(",", ".")) if entry5.get() else None
                C = float(entry6.get().replace(",", ".")) if entry6.get() else None
                
                A, B, C, a, b, c = lei_dos_cossenos(a, b, c, A, B, C)
                resultado.set(f"Ângulo A = {A:.2f}°\nÂngulo B = {B:.2f}°\nÂngulo C = {C:.2f}°\nLado a = {a:.4f}\nLado b = {b:.4f}\nLado c = {c:.4f}")

        elif opcao == "Álgebra":
            escolha = op_alg.get()
            if escolha == "Progressão Aritmética":
                a1 = float(entry1.get().replace(",", "."))
                n = int(entry2.get())
                r = float(entry3.get().replace(",", "."))
                an, sn = progressao_aritmetica(a1, n, r)
                resultado.set(f"Termo a{n} = {an:.4f}\nSoma dos {n} termos = {sn:.4f}")
            elif escolha == "Progressão Geométrica":
                a1 = float(entry1.get().replace(",", "."))
                n = int(entry2.get())
                q = float(entry3.get().replace(",", "."))
                an, sn = progressao_geometrica(a1, n, q)
                resultado.set(f"Termo a{n} = {an:.4f}\nSoma dos {n} termos = {sn:.4f}")
            elif escolha == "Análise Combinatória":
                n = int(entry1.get())
                k = int(entry2.get()) if entry2.get() else None
                if k is None:
                    res = permutacao(n)
                    resultado.set(f"P({n}) = {res}")
                else:
                    comb = combinacao(n, k)
                    perm = permutacao(n, k)
                    arranj = arranjo(n, k)
                    resultado.set(f"C({n},{k}) = {comb}\nP({n},{k}) = {perm}\nA({n},{k}) = {arranj}")

        elif opcao == "Probabilidade":
            escolha = op_prob.get()
            if escolha == "Probabilidade Simples":
                favoraveis = float(entry1.get().replace(",", "."))
                possiveis = float(entry2.get().replace(",", "."))
                prob = probabilidade_evento(favoraveis, possiveis)
                resultado.set(f"P = {prob:.4f} ou {prob*100:.2f}%")
            elif escolha == "Probabilidade Condicional":
                p_a_e_b = float(entry1.get().replace(",", "."))
                p_b = float(entry2.get().replace(",", "."))
                prob = probabilidade_condicional(p_a_e_b, p_b)
                resultado.set(f"P(A|B) = {prob:.4f} ou {prob*100:.2f}%")
            elif escolha == "Teorema de Bayes":
                p_b_dado_a = float(entry1.get().replace(",", "."))
                p_a = float(entry2.get().replace(",", "."))
                p_b = float(entry3.get().replace(",", "."))
                prob = probabilidade_bayes(p_b_dado_a, p_a, p_b)
                resultado.set(f"P(A|B) = {prob:.4f} ou {prob*100:.2f}%")

        elif opcao == "Matrizes":
            escolha = op_matriz.get()
            
            def parse_matriz(texto):
                linhas = texto.split(';')
                matriz = []
                for linha in linhas:
                    elementos = [float(x.replace(",", ".")) for x in linha.split()]
                    matriz.append(elementos)
                return np.array(matriz) if NUMPY_DISPONIVEL else matriz
            
            if escolha in ["Soma", "Subtração", "Multiplicação"]:
                matriz_a = parse_matriz(entry1.get())
                matriz_b = parse_matriz(entry2.get())
                
                if escolha == "Soma":
                    resultado_matriz = matriz_soma(matriz_a, matriz_b)
                elif escolha == "Subtração":
                    resultado_matriz = matriz_subtracao(matriz_a, matriz_b)
                elif escolha == "Multiplicação":
                    resultado_matriz = matriz_multiplicacao(matriz_a, matriz_b)
                
                resultado.set(f"Resultado:\n{resultado_matriz}")
                
            elif escolha == "Transposta":
                matriz = parse_matriz(entry1.get())
                resultado_matriz = matriz_transposta(matriz)
                resultado.set(f"Matriz transposta:\n{resultado_matriz}")
                
            elif escolha == "Determinante":
                matriz = parse_matriz(entry1.get())
                det = matriz_determinante(matriz)
                resultado.set(f"Determinante = {det:.4f}")
                
            elif escolha == "Inversa":
                matriz = parse_matriz(entry1.get())
                inversa = matriz_inversa(matriz)
                if inversa is None:
                    resultado.set("Matriz singular (não possui inversa)")
                else:
                    resultado.set(f"Matriz inversa:\n{inversa}")

        elif opcao == "Sistemas Lineares":
            if not NUMPY_DISPONIVEL:
                resultado.set("Funcionalidade requer NumPy. Instale com: pip install numpy")
            else:
                n = int(entry1.get())
                resultado.set(f"Para {n} variáveis, digite a matriz de coeficientes\n e o vetor de constantes separados por ponto e vírgula")

        elif opcao == "Estatística":
            dados = entry1.get()
            media, mediana, moda, variancia, desvio_padrao = medidas_estatisticas(dados)
            resultado.set(f"Média = {media:.4f}\nMediana = {mediana:.4f}\nModa = {moda}\nVariância = {variancia:.4f}\nDesvio Padrão = {desvio_padrao:.4f}")

        elif opcao == "Matemática Financeira":
            escolha = op_fin.get()
            if escolha == "Juros Simples":
                C = float(entry1.get().replace(",", "."))
                i = float(entry2.get().replace(",", "."))
                n = float(entry3.get().replace(",", "."))
                j, m = juros_simples(C, i, n)
                resultado.set(f"Juros = R$ {j:.2f}\nMontante = R$ {m:.2f}")
            elif escolha == "Prestação Fixa":
                valor = float(entry1.get().replace(",", "."))
                taxa = float(entry2.get().replace(",", "."))
                periodos = int(entry3.get())
                prestacao = valor_prestacao_parcela(valor, taxa, periodos)
                resultado.set(f"Prestação = R$ {prestacao:.2f}\nTotal pago = R$ {prestacao * periodos:.2f}")

        elif opcao == "Cálculo":
            coeficientes = [float(x) for x in entry1.get().replace(",", ".").split()]
            derivada = derivada_polinomial(coeficientes)
            integral = integral_polinomial(coeficientes)
            
            # Formatando os resultados
            grau_original = len(coeficientes) - 1
            grau_derivada = len(derivada) - 1
            grau_integral = len(integral) - 1
            
            func_original = " + ".join([f"{coeficientes[i]}x^{grau_original-i}" for i in range(len(coeficientes))])
            func_derivada = " + ".join([f"{derivada[i]}x^{grau_derivada-i}" for i in range(len(derivada))])
            func_integral = " + ".join([f"{integral[i]}x^{grau_integral-i}" for i in range(len(integral))])
            
            resultado.set(f"Função original: {func_original}\nDerivada: {func_derivada}\nIntegral: {func_integral} + C")

        elif opcao == "Calculadora normal":
            expr = entry1.get()
            expr = expr.replace(",", ".").replace("^", "**")
            res = eval(expr)
            resultado.set(f"Resultado = {res}")

    except ValueError:
        messagebox.showerror("Erro", "Digite valores numéricos válidos!")
    except ZeroDivisionError:
        messagebox.showerror("Erro", "Divisão por zero não é permitida!")
    except Exception as e:
        messagebox.showerror("Erro", f"Ocorreu um erro:\n{str(e)}")

# --------- Interface Tkinter ---------
janela = tk.Tk()
janela.title("Calculadora Matemática Completa")
janela.geometry("700x600")  # Reduzida a altura para caber melhor na tela
janela.resizable(False, False)
janela.configure(bg="#2c3e50")

# Desabilita o foco automático em campos de entrada
janela.option_add('*Entry*TakeFocus', 0)

# Frame principal com scrollbar
main_frame = tk.Frame(janela, bg="#2c3e50")
main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)  # Reduzido o padding vertical

# Canvas e Scrollbar para permitir rolagem
canvas = tk.Canvas(main_frame, bg="#2c3e50", height=500)  # Altura fixa para o canvas
scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
scrollable_frame = tk.Frame(canvas, bg="#2c3e50")

scrollable_frame.bind(
    "<Configure>",
    lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
)

canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
canvas.configure(yscrollcommand=scrollbar.set)

canvas.pack(side="left", fill="both", expand=True)
scrollbar.pack(side="right", fill="y")

# Título
tk.Label(scrollable_frame, text="Calculadora Matemática Completa", font=("Arial", 14, "bold"), 
         bg="#2c3e50", fg="#ecf0f1").pack(pady=3)  # Reduzido o padding

# Área de resultado (em cima) - REDUZIDA
frame_resultado = tk.Frame(scrollable_frame, bg="#34495e", relief=tk.RAISED, bd=1)  # Borda mais fina
frame_resultado.pack(fill=tk.X, pady=3)  # Menor padding

resultado = tk.StringVar()
lbl_resultado = tk.Label(frame_resultado, textvariable=resultado, justify="left", 
                        font=("Consolas", 8), bg="#34495e", fg="#ecf0f1", 
                        height=2, wraplength=650)  # Reduzida a altura e fonte
lbl_resultado.pack(padx=8, pady=3)  # Menor padding

# Preview para calculadora normal
preview = tk.StringVar()
lbl_preview = tk.Label(frame_resultado, textvariable=preview, justify="right",
                      font=("Consolas", 8), bg="#34495e", fg="#bdc3c7",
                      height=1)
lbl_preview.pack(anchor="e", padx=8, pady=(0, 3))  # Menor padding

# Frame para operações
frame_operacoes = tk.Frame(scrollable_frame, bg="#2c3e50")
frame_operacoes.pack(fill=tk.X, pady=3)  # Menor padding

# Escolha da operação
tk.Label(frame_operacoes, text="Escolha a operação:", bg="#2c3e50", 
         fg="#ecf0f1", font=("Arial", 11)).pack()  # Fonte ligeiramente menor

combo = ttk.Combobox(frame_operacoes, values=[
    "Equação 2º Grau",
    "Círculo",
    "Pitágoras",
    "Juros compostos",
    "Polígonos Regulares",
    "Prismas/Sólidos",
    "Trigonometria",
    "Álgebra",
    "Probabilidade",
    "Matrizes",
    "Sistemas Lineares",
    "Estatística",
    "Matemática Financeira",
    "Cálculo",
    "Calculadora normal"
], state="readonly", width=22, font=("Arial", 9))  # Menor largura e fonte
combo.current(0)
combo.pack(pady=3)  # Menor padding
combo.bind("<<ComboboxSelected>>", atualizar_labels)

# Frame para entradas
frame_entradas = tk.Frame(frame_operacoes, bg="#2c3e50")
frame_entradas.pack(pady=3)  # Menor padding

# Widgets principais
label1 = tk.Label(frame_entradas, text="a:", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))  # Fonte menor
label1.pack()
entry1 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))  # Reduzido o tamanho e fonte
entry1.pack(pady=1)  # Menor padding
entry1.bind("<KeyRelease>", atualizar_previa)
# Desabilita o foco automático e o teclado virtual
entry1.configure(takefocus=0)

label2 = tk.Label(frame_entradas, text="b:", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))
entry2 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))
entry2.configure(takefocus=0)

label3 = tk.Label(frame_entradas, text="c:", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))
entry3 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))
entry3.configure(takefocus=0)

label4 = tk.Label(frame_entradas, text="d:", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))
entry4 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))
entry4.configure(takefocus=0)

label5 = tk.Label(frame_entradas, text="e:", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))
entry5 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))
entry5.configure(takefocus=0)

label6 = tk.Label(frame_entradas, text="f:", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))
entry6 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))
entry6.configure(takefocus=0)

# Combos internos para operações específicas
op_circ = ttk.Combobox(frame_entradas, values=["Área", "Circunferência", "Diámetro"], 
                      state="readonly", width=16, font=("Arial", 8))  # Reduzido o tamanho
op_circ.current(0)

op_pyth = ttk.Combobox(frame_entradas, values=["Hipotenusa", "Cateto faltante"], 
                      state="readonly", width=16, font=("Arial", 8))
op_pyth.current(0)
op_pyth.bind("<<ComboboxSelected>>", lambda e: atualizar_pitagoras())

op_poligono = ttk.Combobox(frame_entradas, values=["Área", "Ângulos"], 
                          state="readonly", width=16, font=("Arial", 8))
op_poligono.current(0)
op_poligono.bind("<<ComboboxSelected>>", lambda e: atualizar_poligono())

op_prisma = ttk.Combobox(frame_entradas, values=["Prisma Retangular", "Prisma Triangular", "Cilindro", "Esfera", "Cone"], 
                        state="readonly", width=16, font=("Arial", 8))
op_prisma.current(0)
op_prisma.bind("<<ComboboxSelected>>", lambda e: atualizar_prisma())

op_trig = ttk.Combobox(frame_entradas, values=["Lei dos Senos", "Lei dos Cossenos"], 
                      state="readonly", width=16, font=("Arial", 8))
op_trig.current(0)
op_trig.bind("<<ComboboxSelected>>", lambda e: atualizar_trigonometria())

op_alg = ttk.Combobox(frame_entradas, values=["Progressão Aritmética", "Progressão Geométrica", "Análise Combinatória"], 
                      state="readonly", width=16, font=("Arial", 8))
op_alg.current(0)
op_alg.bind("<<ComboboxSelected>>", lambda e: atualizar_algebra())

op_prob = ttk.Combobox(frame_entradas, values=["Probabilidade Simples", "Probabilidade Condicional", "Teorema de Bayes"], 
                      state="readonly", width=16, font=("Arial", 8))
op_prob.current(0)
op_prob.bind("<<ComboboxSelected>>", lambda e: atualizar_probabilidade())

op_matriz = ttk.Combobox(frame_entradas, values=["Soma", "Subtração", "Multiplicação", "Transposta", "Determinante", "Inversa"], 
                      state="readonly", width=16, font=("Arial", 8))
op_matriz.current(0)
op_matriz.bind("<<ComboboxSelected>>", lambda e: atualizar_matrizes())

op_fin = ttk.Combobox(frame_entradas, values=["Juros Simples", "Prestação Fixa"], 
                     state="readonly", width=16, font=("Arial", 8))
op_fin.current(0)
op_fin.bind("<<ComboboxSelected>>", lambda e: atualizar_financeira())

# Campos específicos para Pitágoras
label_pyth1 = tk.Label(frame_entradas, text="", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))
entry_pyth1 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))
entry_pyth1.configure(takefocus=0)

label_pyth2 = tk.Label(frame_entradas, text="", bg="#2c3e50", fg="#ecf0f1", font=("Arial", 9))
entry_pyth2 = tk.Entry(frame_entradas, width=22, font=("Arial", 9))
entry_pyth2.configure(takefocus=0)

# Botões de ação
frame_botoes = tk.Frame(frame_operacoes, bg="#2c3e50")
frame_botoes.pack(pady=3)  # Menor padding

btn_calcular = tk.Button(frame_botoes, text="Calcular", font=("Arial", 9, "bold"), 
                bg="#27ae60", fg="white", width=8, command=calcular)  # Botão menor
btn_calcular.pack(side=tk.LEFT, padx=3)  # Menor padding

btn_limpar = tk.Button(frame_botoes, text="Limpar Tudo", font=("Arial", 9), 
                bg="#e74c3c", fg="white", width=8, command=limpar_tudo)  # Botão menor
btn_limpar.pack(side=tk.LEFT, padx=3)  # Menor padding

# --------- Teclado Numérico Ampliado ---------
frame_teclado = tk.Frame(scrollable_frame, bg="#2c3e50")
frame_teclado.pack(pady=3, anchor="w")  # Menor padding e alinhado à esquerda

tk.Label(frame_teclado, text="Teclado:", font=("Arial", 11, "bold"), 
         bg="#2c3e50", fg="#ecf0f1").pack()

# Frame para teclas principais
frame_teclas_principais = tk.Frame(frame_teclado, bg="#2c3e50")
frame_teclas_principais.pack()

# Teclas numéricas (botões menores)
botoes_numericos = [
    ['7', '8', '9', '(', ')', '⌫'],
    ['4', '5', '6', '*', '/', 'C'],
    ['1', '2', '3', '+', '-', '='],
    ['0', '.', ',', '^', 'π', '√']
]

for linha in botoes_numericos:
    frame_linha = tk.Frame(frame_teclas_principais, bg="#2c3e50")
    frame_linha.pack()
    for tecla in linha:
        if tecla == '=':
            btn = tk.Button(frame_linha, text=tecla, width=3, height=1,  # Botões menores
                           font=("Arial", 8, "bold"), bg="#3498db", fg="white",
                           command=calcular)
        elif tecla == 'C':
            btn = tk.Button(frame_linha, text=tecla, width=3, height=1,
                           font=("Arial", 8, "bold"), bg="#e74c3c", fg="white",
                           command=limpar_tudo)
        elif tecla == '⌫':
            btn = tk.Button(frame_linha, text=tecla, width=3, height=1,
                           font=("Arial", 8), bg="#95a5a6", fg="white",
                           command=apagar_caractere)
        elif tecla in ['π', '√']:
            btn = tk.Button(frame_linha, text=tecla, width=3, height=1,
                           font=("Arial", 8, "bold"), bg="#16a085", fg="white",
                           command=lambda t=tecla: inserir_texto('math.pi' if t == 'π' else 'math.sqrt('))
        else:
            btn = tk.Button(frame_linha, text=tecla, width=3, height=1,
                           font=("Arial", 8), bg="#34495e", fg="#ecf0f1",
                           command=lambda t=tecla: inserir_texto(t))
        btn.pack(side=tk.LEFT, padx=1, pady=1)  # Menor espaçamento

# Teclas de funções especiais (botões menores)
frame_funcoes = tk.Frame(frame_teclado, bg="#2c3e50")
frame_funcoes.pack(pady=2, anchor="w")  # Menor padding e alinhado à esquerda

funcoes_especiais = [
    ('sin', 'math.sin('),
    ('cos', 'math.cos('),
    ('tan', 'math.tan('),
    ('log', 'math.log10('),
    ('ln', 'math.log('),
    ('°', 'math.radians('),
    ('rad', 'math.degrees(')
]

for i in range(0, len(funcoes_especiais), 4):
    frame_linha = tk.Frame(frame_funcoes, bg="#2c3e50")
    frame_linha.pack()
    for funcao in funcoes_especiais[i:i+4]:
        texto, comando = funcao
        btn = tk.Button(frame_linha, text=texto, width=5, height=1,  # Botões menores
                       font=("Arial", 7), bg="#8e44ad", fg="white",
                       command=lambda c=comando: inserir_texto(c))
        btn.pack(side=tk.LEFT, padx=1, pady=1)  # Menor espaçamento

# Inicializa labels corretos
atualizar_labels()

janela.mainloop()
