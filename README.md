Feito por Álvaro Carvalho M. de Almeida

-------------------------------------------------------------------------------------------------------------------------------------------------------------------

Resolução numérica da equação de Laplace na forma discretizada, como forma de simular o processo de condução de calor em regime permanente 2D.
Além disso, foi utilizado o método iterativo de Liebmann/Gauss-Seidel, que se mostrou mais eficiente no quesito custo computacional (malhas muito finas representariam matrizes imensas).
Todo o processo foi norteado pela literatura de Chapra & Canale (Métodos Numéricos para Engenharia, 5ed.) e pelos vídeos da discretização passo a passo do Prof. Rafael Gabler Gontijo.
Há também a análise do fluxo de calor resultante de materiais como o Alumínio e Níquel puros, através do coeficiente de condutividade térmica (k) pela Lei de Fourier de condução.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------

EQUAÇÃO DE LAPLACE PARA CONDUÇÃO DE CALOR 2D

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

OBS:REGIME PERMANENTE

FORMA DISCRETIZADA DA EQUAÇÃO

$$T_{i,j} = \frac{T_{i+1,j} + T_{i-1,j} + T_{i,j+1} + T_{i,j-1}}{4}$$

LEI DE FOURIER PARA O FLUXO DE CALOR

$$\vec{q}'' = -k \nabla T$$

OBS:REGIME PERMANENTE = GRADIENTE DE TEMPERATURA IGUAL

-------------------------------------------------------------------------------------------------------------------------------------------------------------------

Observação:
Há alternativas open-source do MATLAB, como o Octave (que aceitam arquivos .m).

-------------------------------------------------------------------------------------------------------------------------------------------------------------------
