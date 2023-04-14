# Otimizacao-Topologica

A Otimização Topológica consiste num método computacional poderoso de otimização estrutural que permite projetar a topologia ótima de estruturas segundo um certo critério de custo. O algoritmo do Método de Otimização Topológica (MOT) combina métodos de otimização com um método numérico de análise, como o Método dos Elementos Finitos (MEF), para distribuir material no interior de um domínio de projeto de maneira a minimizar (ou maximizar) a função objetivo especificada, satisfazendo dadas restrições impostas ao problema de otimização. O principal objetivo deste projeto de iniciação científica é explorar a utilização da linguagem Python na implementação computacional do algoritmo do MOT utilizando o modelo SIMP (“Solid Isotropic Material with Penalization”). Pretende-se desenvolver uma ferramenta computacional (em Python) de baixo custo capaz de lidar com modelos 3D para projeto de estruturas aeroespaciais. No projeto, o problema será estudado através de implementação de rotinas computacionais de MEF e de otimização topológica, bem como modelagem e simulação computacional das estruturas obtidas.

# Alguns resultados do programa.

![image](https://user-images.githubusercontent.com/128917882/232026190-1e9aa926-8557-45d8-84c1-d48fcd25ffae.png) ![image](https://user-images.githubusercontent.com/128917882/232026209-bf23df3a-c9b8-45d5-9d44-aaa1cfce7856.png)

![image](https://user-images.githubusercontent.com/128917882/232026256-51025959-c029-4f0c-a631-b19574f118ce.png) ![image](https://user-images.githubusercontent.com/128917882/232026290-00f5a130-f577-4cd4-8aa7-cfd5cd2cd22c.png)

![image](https://user-images.githubusercontent.com/128917882/232026321-16fe3b6c-8e01-42e0-9543-37b91b6850c9.png) ![image](https://user-images.githubusercontent.com/128917882/232026459-dbf547bf-4b57-460b-97fa-b73bf43caa06.png)

# Explicação do código

Nas linhas 1 a 9, bibliotecas instaladas previamente são importadas para o script para que se possa utilizar suas funções. A biblioteca “time” fornece meios de calcular o tempo que o programa demora para ir de uma linha inicial até uma linha final.

![image](https://user-images.githubusercontent.com/128917882/232022985-5f8af64b-cd4c-4074-b727-1d87add2f0b7.png)

Nas linhas 11 a 24, inicia-se o programa principal através de uma função que posteriormente poderá ser chamada pelo usuário para realizar a otimização da estrutura. É dentro dessa função, todos os próximos segmentos do código estarão inseridos. Além disso, tem-se dois parâmetros do looping que será estabelecido posteriormente: “maxloop”, que é o valor máximo de iterações que o programa pode atingir; “tolx”, que é o valor do critério de parada do programa. Nesse trecho também são inseridas as propriedades intrínsecas do material, como já descrito na seção anterior.

![image](https://user-images.githubusercontent.com/128917882/232023113-20f44e71-0cab-4664-ab26-98f6461bdeec.png)

Para entender as próximas partes da implementação é utilizado as notações e convenções feitas por Tovar et al. 2014 [14]. Primeiramente, é necessário enumerar os nós dos elementos da malha, que será chamado de ID do nó (NID). Essa numeração será feita de cima para baixo e de trás para frente, assim como mostra o exemplo da figura 22.

![image](https://user-images.githubusercontent.com/128917882/232024676-1da1ca07-c470-40a7-b2f8-92217c489f6f.png)

Separando um elemento qualquer dessa estrutura, tem-se um cubo de lado unitário com coordenadas x, y e z, como mostrado na figura 23.

![image](https://user-images.githubusercontent.com/128917882/232024756-80ed0b1e-9b24-4cba-a2ee-3ad667417e76.png)

Cada nó da figura 22, tem 3 graus de liberdade e por convenção a identificação de cada grau de liberdade é feita através, também, de uma numeração. Essa numeração segue a seguinte sequência: NID1 (Corresponde ao nó 1 da figura 22) tem os graus de liberdade 1 na direção x, 2 na direção y e 3 na direção z; NID2 tem os graus de liberdade 4 na direção x, 5 na direção y e 6 na direção z; e assim por diante.

As convenções acima podem ser generalizadas e resumidas na tabela 1. Com isso é possível implementar as linhas 25 a 60 do programa.

![image](https://user-images.githubusercontent.com/128917882/232024838-272618e4-e72b-4cdf-a47e-d03b5a59552f.png)

![image](https://user-images.githubusercontent.com/128917882/232024980-2255f33f-1f8a-4c45-96a5-94f15bd4f76e.png)

![image](https://user-images.githubusercontent.com/128917882/232025068-1159d590-2be6-4a6a-b947-091c44dac925.png)


No trecho acima, constrói-se a malha da estrutura. Primeiro define-se o local de aplicação da força, através dos Ids dos nós, e a direção que a força através da convenção de numeração dos graus de liberdade. Posteriormente, define-se os nós e graus de liberdade que estão fixos. E assim, tem-se o domínio ꭥ. 

Nesse trecho também, é calculado a matriz de rigidez de um elemento “KE” chamando a função “lk_H8” que será explicada mais a frente. Por fim, a matriz “edofMat” e os valores “iK” e “jK”, que já foram explicados na seção anterior, foram calculados. 

As linhas 61 a 88 são equivalentes as linhas 46 a 64 do programa para estruturas 2D da seção anterior. As únicas mudanças realizadas são devidas ao aumento de uma direção no problema, assim é necessário adicionar mais um laço “for” e a distância Δ(e,i) tem uma nova componente.

![image](https://user-images.githubusercontent.com/128917882/232025152-b8b62a86-5103-4bc7-8639-7329a9a218b3.png)

![image](https://user-images.githubusercontent.com/128917882/232025242-e804fece-aa78-468d-8661-44d65acc7e02.png)

No caso das próximas linhas, 90 a 95, são apenas inicializações de variáveis que serão utilizadas posteriormente, nesse caso são vetores com apenas o número 1 em todas as posições.

![image](https://user-images.githubusercontent.com/128917882/232025294-cb11c9cc-9dd6-4c12-a2f1-2e9a4a3783a5.png)

As linhas 97 a 140, tem os mesmos conceitos das linhas 80 a 113 do código para estruturas 2D, exceto a solução do sistema. Para resolver problemas em 3D necessitamos de quantidades elevadas de elementos, isso gera problemas computacionais para resolução da matriz inversa na restrição de equilíbrio na equação (1), tais como: necessidade de muita memória RAM para armazenar dados e a velocidade de processamento fica prejudicada.
Para contornar esses problemas, foi utilizado o método dos gradientes conjugados (método cg). O método cg é um método iterativo que acaba quando atinge uma tolerância de erro (“tol”) ou em um número máximo de iterações. Iniciando com um valor x_0, novas estimativas x_1, x_2, e assim por diante, são determinadas, cada vez mais próximas da solução real. Mais detalhes do método são encontrados em [15].

![image](https://user-images.githubusercontent.com/128917882/232025544-89ebd404-6cad-458d-a894-108c926f96f2.png)

![image](https://user-images.githubusercontent.com/128917882/232025616-eb78c55b-e88f-497b-9879-1af201e644e2.png)

![image](https://user-images.githubusercontent.com/128917882/232025642-3a64d7b1-d62d-4adc-8259-4abf4f4ef5b4.png)

No segmento a seguir é plotado a figura em 3D, para isso foi utilizado o a estrutura de “voxels” existente no Python. A estratégia aqui é passar por todo o vetor de densidade “xPhys”, lembrando das convenções anteriores, para plotar esses cubos sempre que a densidade do elemento for maior que 0.5 e traçando uma escala de cinza através da equação da linha 160.

![image](https://user-images.githubusercontent.com/128917882/232025721-0cbb36e1-4218-44cb-a62a-fc7333178683.png)

Nas linhas 165 a 216, uma função é construída com o intuito de gerar a matriz de rigidez de um elemento, ou seja, essa função retorna k_e^0.

![image](https://user-images.githubusercontent.com/128917882/232025840-89280066-b235-4f1b-a895-d105eb581d4e.png)

![image](https://user-images.githubusercontent.com/128917882/232025877-6cdcd12f-c81e-496b-af27-7909d27ebbaf.png)

![image](https://user-images.githubusercontent.com/128917882/232025906-4e1c38e4-2959-4bb5-b564-5328693ee04f.png)

Nas linhas 218 a 233, é calculado o novo valor de densidade (xnew) através do critério de optimalidade dado pela equação (2), assim como mostrado no programa para estruturas bidimensionais. Porém, para estruturas tridimensionais, o vetor de densidade tem “nelx*nely*nelz” elementos.

![image](https://user-images.githubusercontent.com/128917882/232025977-d40b2605-1b15-47f8-8338-431d9bd46276.png)

Por fim, as últimas linhas abaixo apenas funcionam se não digitarmos no terminal a função de chamada “topopt3d(nelx,nely,nelz,volfrac,penal,rmin,ft)”. Assim podemos mudar os valores da variável direto no código e rodar diretamente na IDE.

![image](https://user-images.githubusercontent.com/128917882/232026078-4b86d5b0-3718-4777-a7f4-9ae572e5ba35.png)
