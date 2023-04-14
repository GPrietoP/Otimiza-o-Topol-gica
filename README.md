# Otimizacao-Topologica

A Otimização Topológica consiste num método computacional poderoso de otimização estrutural que permite projetar a topologia ótima de estruturas segundo um certo critério de custo. O algoritmo do Método de Otimização Topológica (MOT) combina métodos de otimização com um método numérico de análise, como o Método dos Elementos Finitos (MEF), para distribuir material no interior de um domínio de projeto de maneira a minimizar (ou maximizar) a função objetivo especificada, satisfazendo dadas restrições impostas ao problema de otimização. O principal objetivo deste projeto de iniciação científica é explorar a utilização da linguagem Python na implementação computacional do algoritmo do MOT utilizando o modelo SIMP (“Solid Isotropic Material with Penalization”). Pretende-se desenvolver uma ferramenta computacional (em Python) de baixo custo capaz de lidar com modelos 3D para projeto de estruturas aeroespaciais. No projeto, o problema será estudado através de implementação de rotinas computacionais de MEF e de otimização topológica, bem como modelagem e simulação computacional das estruturas obtidas.

Nas linhas 1 a 9, bibliotecas instaladas previamente são importadas para o script para que se possa utilizar suas funções. A biblioteca “time” fornece meios de calcular o tempo que o programa demora para ir de uma linha inicial até uma linha final.

![image](https://user-images.githubusercontent.com/128917882/232022985-5f8af64b-cd4c-4074-b727-1d87add2f0b7.png)

Nas linhas 11 a 24, inicia-se o programa principal através de uma função que posteriormente poderá ser chamada pelo usuário para realizar a otimização da estrutura. É dentro dessa função, todos os próximos segmentos do código estarão inseridos. Além disso, tem-se dois parâmetros do looping que será estabelecido posteriormente: “maxloop”, que é o valor máximo de iterações que o programa pode atingir; “tolx”, que é o valor do critério de parada do programa. Nesse trecho também são inseridas as propriedades intrínsecas do material, como já descrito na seção anterior.

![image](https://user-images.githubusercontent.com/128917882/232023113-20f44e71-0cab-4664-ab26-98f6461bdeec.png)





