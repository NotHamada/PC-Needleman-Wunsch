/*

Programa de Alinhamento de Genoma Versao 27/07/2024.

Programa "sequencial" baseado no metodo Needleman-Wunsch, desenvolvido por Saul
B. Needleman e Christian D. Wunsch, publicado em 1970, conforme segue:

Needleman, Saul B. & Wunsch, Christian D. (1970)."A general method applicable to
the search for similarities in the amino acid sequence of two proteins". Journal
of Molecular Biology.48 (3):443-53.

Este programa NAO pode ser usado, da forma original ou modificada, completa ou
parcial, para finalidade comercial ou lucrativa, sem a previa autorizacao dos
seus desenvolvedores, os quais detem os direitos autorais.

Este programa PODE ser livremente e gratuitamente usado, da forma original ou
modificada, completa ou parcial, como apoio ao processo de ensino-aprendizagem,
nao comercial e nao lucrativo, por qualquer pessoa, desde que os resultados
obtidos ou os produtos/subprodutos gerados tambem possam ser usados com a mesma
liberdade e gratuidade.

Todos os desenvolvedores responsaveis pelo codigo original e futuras modificacoes
devem ser informados na lista a seguir, apos o ultimo, antes da distribuicao do
codigo modificado, informando quais foram as modificacoes realizadas.

Lista de Desenvolvedores:

Desenvolvedor: Ronaldo Augusto de Lara Goncalves
Data da ultima atualizacao: 27/07/2024
eMail: ralgonca@gmail.com
Whatsapp: 55(44)99159-2310
Modificacoes realizadas: desenvolvimento do codigo original, composto pelos
modulos leTamMaior(), leTamMenor(), lePenalidade(), menuOpcao(void), trataOpcao(),
geraSequencias(), leMatrizPesos(), mostraMatrizPesos(), leGrauMutacao(),
geraMatrizEscores(), mostraMatrizEscores(), leSequencias(), geraSequencias(),
mostraSequencias(), traceBack(), mostraAlinhamentoGlobal() e main().



======================================
BREVE DESCRICAO:

======================================

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define A 0 // representa uma base Adenina
#define T 1 // representa uma base Timina
#define G 2 // representa uma base Guanina
#define C 3 // representa uma base Citosina
#define X 4 // representa um gap

#define sair 11

/* mapaBases mapeia indices em caracteres que representam as bases, sendo 0='A',
1='T', 2='G', 3='C' e 4='-' representando gap */

char mapaBases[5] = {'A', 'T', 'G', 'C', '-'};

/* seqMaior e seqMenor representam as duas sequencias de bases de entrada, a
   serem comparadas, inicializadas conforme segue. Elas conterao os indices aos
   inves dos proprios caracteres. seqMenor deve ser menor ou igual a seqMaior. */

int maxSeq = 1000; // tamanho maximo de bases em uma sequencia genomica

int seqMaior[1000] = {A, A, C, T, T, A},
    seqMenor[1000] = {A, C, T, T, G, A};

/* alinhaGMaior representa a sequencia maior ja alinhada, assim como alinhaGMenor,
   ambas obtidas no traceback. As duas juntas, pareadas, formam o alinhamento
   global. Tal alinhamento global pode ser obtido de duas formas: a partir do
   primeiro maior escore ou a partir do ultimo maior escore */

int alinhaGMaior[1000],
    alinhaGMenor[1000];

/* matrizEscores representa a matriz de escores que sera preenchida pelo metodo.
   A matriz, ao final de seu preenchimento, permitira obter o melhor alinhamento
   global entre as sequencias seqMaior e seqMenor, por meio de uma operacao
   denominada TraceBack. Uma linha e uma coluna extras sao adicionadas na matriz
   para inicializar as pontuacoes/escores. Trata-se da linha 0 e coluna 0. A
   matriz de escores tera tamSeqMenor+1 linhas e tamSeqMaior+1 colunas.
   Considera-se a primeira dimensao da matriz como linhas e a segunda como colunas.*/

int matrizEscores[1000 + 1][1000 + 1];

int tamSeqMaior = 6, /* tamanho da sequencia maior, inicializado como 6 */
    tamSeqMenor = 6, /* tamanho da sequencia menor, inicializado como 6 */
    tamAlinha,       /* tamanho do alinhamento global obtido */
    penalGap = 0,    /* penalidade de gap, a ser descontada no escore acumulado
                        quando um gap eh encontrado */
    grauMuta = 0,    /* porcentagem maxima de mutacao na geracao aleatoria da
                        sequencia menor, a qual eh copiada da maior e sofre algumas
                        trocas de bases */
    escoreDiag,      /* escore da diagonal anterior da matriz de escores */
    escoreLin,       /* escore da linha anterior da matriz de escores */
    escoreCol,
    prob,
    resp_geracao; /* escore da coluna anterior da matriz de escores */

/*  matrizPesos contem os pesos do pareamento de bases. Estruturada e inicializada
    conforme segue, onde cada linha ou coluna se refere a uma das bases A, T, G
    ou C. Considera-se a primeira dimensao da matriz como linhas e a segunda como
    colunas. Na configuracao default, o peso de bases iguais eh 1 e o peso de bases
    diferentes eh 0, mas pode ser alterado. Outra configuracao usual eh 2 para
    bases iguais e -1 para bases diferentes

0 1 2 3 A T G C 0 A 1 0 0 0 1 T 0 1 0 0 2 G 0 0 1 0 3 C 0 0 0 1 */

int matrizPesos[4][4] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

int indRef = -1,                  // indice da sequencia maior a partir do qual extrai a sequencia
                                  // menor, no caso de geracao aleatoria
    nTrocas = -1,                 // quantidade de trocas na geracao automatica da sequencia menor,
                                  // a partir de um segmento da sequencia maior
    linPMaior, colPMaior, PMaior, // suporte para deteccao do primeiro maior escore
    linUMaior, colUMaior, UMaior; // suporte para deteccao do ultimo maior escore

/* leitura do tamanho da sequencia maior */

void leTamMaior(int rank)
{
  if (rank == 0)
  {
    printf("\nLeitura do Tamanho da Sequencia Maior:");
    do
    {
      printf("\nDigite 0 < valor < %d = ", maxSeq);
      scanf("%d", &tamSeqMaior);
    } while ((tamSeqMaior < 1) || (tamSeqMaior > maxSeq));
  }
}

void leTamMenor(int rank)
{
  if (rank == 0)
  {
    printf("\nLeitura do Tamanho da Sequencia Menor:");
    do
    {
      printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
      scanf("%d", &tamSeqMenor);
    } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
  }
}

int lePenalidade(int rank)
{
  int penal;

  if (rank == 0)
  {
    printf("\nLeitura da Penalidade de Gap:");
    do
    {
      printf("\nDigite valor >= 0 = ");
      scanf("%d", &penal);
    } while (penal < 0);

    // Broadcast do valor de penal para todos os processos
    MPI_Bcast(&penal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else
  {
    // Recebe o broadcast do valor de penal
    MPI_Bcast(&penal, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }

  return penal;
}
/* leitura da matriz de pesos */
void leMatrizPesos(int rank)
{
  int i, j;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    // Processo 0 faz a leitura dos dados
    printf("\nLeitura da Matriz de Pesos:\n");
    for (i = 0; i < 4; i++)
    {
      for (j = 0; j < 4; j++)
      {
        printf("Digite valor %c x %c = ", mapaBases[i], mapaBases[j]);
        scanf("%d", &(matrizPesos[i][j]));
      }
      printf("\n");
    }
  }

  // Broadcast da matriz para todos os processos
  MPI_Bcast(matrizPesos, 16, MPI_INT, 0, MPI_COMM_WORLD);
}
/* mostra da matriz de pesos */
void mostraMatrizPesos(void)
{
  int i, j;

  printf("\nMatriz de Pesos Atual:");
  printf("\n%4c%4c%4c%4c%4c\n", ' ', 'A', 'T', 'G', 'C');
  for (i = 0; i < 4; i++)
  {
    printf("%4c", mapaBases[i]);
    for (j = 0; j < 4; j++)
      printf("%4d", matrizPesos[i][j]);
    printf("\n");
  }
}

/* leitura da porcentagem maxima (grau) de mutacao aleatoria. Essa porcentagem eh
   usada na geracao aleatoria da seqMenor. A seqMenor eh obtida a partir da seqMaior, para se parecer com ela, se diferenciando
   por um certo grau de alteracoes em suas bases, fornecida pelo usuario. Esse
   metodo evita a gera��o aleatoria de sequencias totalmente diferentes. A
   quantidade de trocas realizadas eh no maximo a porcentagem aqui informada. */

void leSequenciasDeArquivo(char *fileName, int rank)
{
  FILE *file;
  char buffer[maxSeq + 2]; // Buffer para leitura incluindo o newline e null terminator

  if (rank == 0)
  {
    file = fopen(fileName, "r");
    if (file == NULL)
    {
      printf("Erro ao abrir o arquivo %s.\n", fileName);
      exit(1);
    }

    // Leitura da sequência maior
    if (fgets(buffer, sizeof(buffer), file) != NULL)
    {
      tamSeqMaior = strlen(buffer);
      if (buffer[tamSeqMaior - 1] == '\n')
      {
        buffer[tamSeqMaior - 1] = '\0';
        tamSeqMaior--;
      }
      for (int i = 0; i < tamSeqMaior; i++)
      {
        switch (buffer[i])
        {
        case 'A':
          seqMaior[i] = A;
          break;
        case 'T':
          seqMaior[i] = T;
          break;
        case 'G':
          seqMaior[i] = G;
          break;
        case 'C':
          seqMaior[i] = C;
          break;
        default:
          printf("Caractere inválido na sequência maior: %c\n", buffer[i]);
          fclose(file);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
    }
    else
    {
      printf("Erro ao ler a sequência maior do arquivo %s.\n", fileName);
      fclose(file);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Leitura da sequência menor
    if (fgets(buffer, sizeof(buffer), file) != NULL)
    {
      tamSeqMenor = strlen(buffer);
      if (buffer[tamSeqMenor - 1] == '\n')
      {
        buffer[tamSeqMenor - 1] = '\0';
        tamSeqMenor--;
      }
      for (int i = 0; i < tamSeqMenor; i++)
      {
        switch (buffer[i])
        {
        case 'A':
          seqMenor[i] = A;
          break;
        case 'T':
          seqMenor[i] = T;
          break;
        case 'G':
          seqMenor[i] = G;
          break;
        case 'C':
          seqMenor[i] = C;
          break;
        default:
          printf("Caractere inválido na sequência menor: %c\n", buffer[i]);
          fclose(file);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      }
    }
    else
    {
      printf("Erro ao ler a sequência menor do arquivo %s.\n", fileName);
      fclose(file);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    fclose(file);
  }

  // Broadcast dos tamanhos das sequências para todos os processos
  MPI_Bcast(&tamSeqMaior, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tamSeqMenor, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Após o tamanho ser conhecido, os processos podem ajustar o buffer para as sequências
  // Fazendo o broadcast das sequências agora com os tamanhos corretos
  MPI_Bcast(seqMaior, tamSeqMaior, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(seqMenor, tamSeqMenor, MPI_INT, 0, MPI_COMM_WORLD);
}

int leGrauMutacao(int rank)
{
  if (rank == 0)
  {
    // Somente o processo 0 faz a leitura
    printf("\nLeitura da Porcentagem Maxima de Mutacao Aleatoria:\n");
    do
    {
      printf("\nDigite 0 <= valor <= 100 = ");
      scanf("%d", &prob);
    } while ((prob < 0) || (prob > 100));
  }

  // Broadcast do valor de grauMuta para todos os processos
  MPI_Bcast(&prob, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return prob;
}

/* leitura manual das sequencias de entrada seqMaior e seqMenor */
void leSequencias(int rank)
{
  int i, erro;
  char seqMaiorAux[maxSeq], seqMenorAux[maxSeq];

  if (rank == 0)
  {
    // Leitura da sequência maior
    printf("\nLeitura das Sequencias:\n");
    do
    {
      printf("\nPara a Sequencia Maior,");
      printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
      do
      {
        printf("\n> ");
        fgets(seqMaiorAux, maxSeq, stdin);
        tamSeqMaior = strlen(seqMaiorAux) - 1; // Remove o newline
        seqMaiorAux[tamSeqMaior] = '\0';       // Remove o newline do final da string
      } while (tamSeqMaior < 1);
      printf("\ntamSeqMaior = %d\n", tamSeqMaior);
      i = 0;
      erro = 0;
      do
      {
        switch (seqMaiorAux[i])
        {
        case 'A':
          seqMaior[i] = A;
          break;
        case 'T':
          seqMaior[i] = T;
          break;
        case 'G':
          seqMaior[i] = G;
          break;
        case 'C':
          seqMaior[i] = C;
          break;
        default:
          erro = 1; // Caractere inválido
        }
        i++;
      } while ((erro == 0) && (i < tamSeqMaior));
    } while (erro == 1);

    // Leitura da sequência menor
    do
    {
      printf("\nPara a Sequencia Menor, ");
      printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
      do
      {
        printf("\n> ");
        fgets(seqMenorAux, maxSeq, stdin);
        tamSeqMenor = strlen(seqMenorAux) - 1; // Remove o newline
        seqMenorAux[tamSeqMenor] = '\0';       // Remove o newline do final da string
      } while ((tamSeqMenor < 1) || (tamSeqMenor > tamSeqMaior));
      printf("\ntamSeqMenor = %d\n", tamSeqMenor);

      i = 0;
      erro = 0;
      do
      {
        switch (seqMenorAux[i])
        {
        case 'A':
          seqMenor[i] = A;
          break;
        case 'T':
          seqMenor[i] = T;
          break;
        case 'G':
          seqMenor[i] = G;
          break;
        case 'C':
          seqMenor[i] = C;
          break;
        default:
          erro = 1; // Caractere inválido
        }
        i++;
      } while ((erro == 0) && (i < tamSeqMenor));
    } while (erro == 1);
  }

  // Broadcast dos tamanhos das sequências
  MPI_Bcast(&tamSeqMaior, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tamSeqMenor, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Broadcast das sequências
  MPI_Bcast(seqMaior, tamSeqMaior, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(seqMenor, tamSeqMenor, MPI_INT, 0, MPI_COMM_WORLD);
}

/* geracao das sequencias aleatorias, conforme tamanho. Gera-se numeros aleatorios
   de 0 a 3 representando as bases 'A', 'T', 'G' e 'C'. Gera-se primeiramente a
   maior sequencia e desta extrai a menor sequencia. A menor sequencia eh obtida
   da maior por meio de algumas trocas de bases (mutacoes), de acordo com o grau
   de mutacao informado. A ideia eh gerar sequencias parecidas, mas com certo grau
   de diferenca. */

void geraSequencias(int rank)
{
  int i, dif, probAux;
  char base;

  // Inicializa o gerador de números aleatórios com base no rank do processo
  srand(time(NULL) + rank); // Isso garante que cada processo tenha uma sequência diferente

  if (rank == 0)
  {
    printf("\nGeracao Aleatoria das Sequencias:\n");

    // Gerando a sequência maior
    for (i = 0; i < tamSeqMaior; i++)
    {
      base = rand() % 4; // Produz valores de 0 a 3
      seqMaior[i] = base;
    }

    dif = tamSeqMaior - tamSeqMenor; // Diferença entre os tamanhos das sequências

    indRef = 0;
    if (dif > 0)
      indRef = rand() % dif; // Produz um índice aleatório para indexar a sequência maior

    // Gerando a sequência menor a partir da maior
    for (i = 0; i < tamSeqMenor; i++)
      seqMenor[i] = seqMaior[indRef + i];

    // Causa mutações aleatórias na sequência menor
    i = 0;
    nTrocas = 0;
    while ((i < tamSeqMenor) && (nTrocas < ((grauMuta * tamSeqMenor) / 100)))
    {
      probAux = rand() % 100 + 1;

      if (probAux <= grauMuta)
      {
        seqMenor[i] = (seqMenor[i] + (rand() % 3) + 1) % 4;
        nTrocas++;
      }
      i++;
    }

    printf("\nSequencias Geradas: Dif = %d, IndRef = %d, NTrocas = %d\n", dif, indRef, nTrocas);
  }

  // Broadcast dos tamanhos e sequências para todos os processos
  MPI_Bcast(&tamSeqMaior, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tamSeqMenor, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(seqMaior, tamSeqMaior, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(seqMenor, tamSeqMenor, MPI_INT, 0, MPI_COMM_WORLD);
}

/* mostra das sequencias seqMaior e seqMenor */
void mostraSequencias(void)
{
  int i;

  printf("\nSequencias Atuais:\n");
  printf("\nSequencia Maior, Tam = %d\n", tamSeqMaior);
  for (i = 0; i < tamSeqMaior; i++)
    printf("%c", mapaBases[seqMaior[i]]);
  printf("\n");

  for (i = 0; i < tamSeqMaior; i++)
    if (i != indRef)
      printf(" ");
    else
      printf("^");
  printf("\nIndice de Referencia = %d\n", indRef);

  printf("\nSequencia Menor, Tam = %d\n", tamSeqMenor);
  for (i = 0; i < tamSeqMenor; i++)
    printf("%c", mapaBases[seqMenor[i]]);
  printf("\n");

  for (i = 0; i < tamSeqMenor; i++)
    if (seqMenor[i] != seqMaior[indRef + i])
      printf("^");
    else
      printf(" ");
  printf("\nQuantidade de trocas = %d\n", nTrocas);
}

/* geraMatrizEscores gera a matriz de escores. A matriz de escores tera
   tamSeqMenor+1 linhas e tamSeqMaior+1 colunas. A linha 0 e a coluna
   0 s�o adicionadas para representar gaps e conter penalidades. As
   demais linhas e colunas sao associadas as bases da seqMenor e da
   SeqMaior, respectivamente. */

void geraMatrizEscores(int rank, int size)
{
  int lin, col, peso;
  int escoreDiag, escoreLin, escoreCol;
  int lin_inicial, lin_final;

  // Cada processo calcula uma linha de cada vez
  lin_inicial = rank;
  lin_final = tamSeqMenor + 1;

  // Inicializa a matriz de penalidades no processo 0
  if (rank == 0)
  {
    for (col = 0; col <= tamSeqMaior; col++)
      matrizEscores[0][col] = col * penalGap; // Penalidades de gaps na primeira linha

    for (lin = 0; lin <= tamSeqMenor; lin++)
      matrizEscores[lin][0] = lin * penalGap; // Penalidades de gaps na primeira coluna
  }

  // Envia as penalidades iniciais da matriz para os outros processos
  if (rank == 0)
  {
    for (int dest = 1; dest < size; dest++)
    {
      for (col = 0; col <= tamSeqMaior; col++)
      {
        MPI_Send(matrizEscores[0], tamSeqMaior + 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
      }
    }
  }
  else
  {
    MPI_Recv(matrizEscores[0], tamSeqMaior + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // Processos começam a calcular suas linhas a partir da segunda linha
  for (lin = lin_inicial; lin <= lin_final; lin += size - 1)
  {
    if (rank != 0)
    {
      // Recebe a linha anterior do processo anterior
      MPI_Recv(matrizEscores[lin - 1], tamSeqMaior + 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Calcula a linha atual
    for (col = 1; col <= tamSeqMaior; col++)
    {
      // Acessa o índice correto das bases das sequências
      int baseSeqMenor = seqMenor[lin - 1];
      int baseSeqMaior = seqMaior[col - 1];

      // Obtenha o peso da matriz de pesos
      peso = matrizPesos[baseSeqMenor][baseSeqMaior];

      // Calcula os escores possíveis (diagonal, em cima, à esquerda)
      escoreDiag = matrizEscores[lin - 1][col - 1] + peso;
      escoreLin = matrizEscores[lin - 1][col] - penalGap;
      escoreCol = matrizEscores[lin][col - 1] - penalGap;

      // Escolhe o maior escore
      matrizEscores[lin][col] = escoreDiag;
      if (escoreLin > matrizEscores[lin][col])
        matrizEscores[lin][col] = escoreLin;
      if (escoreCol > matrizEscores[lin][col])
        matrizEscores[lin][col] = escoreCol;
    }

    // Envia a linha para o próximo processo (se houver)
    if (rank != size - 1)
    {
      MPI_Send(matrizEscores[lin], tamSeqMaior + 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
  }

  // O processo 0 coleta todas as linhas e monta a matriz completa
  if (rank == 0)
  {
    for (int p = 1; p < size; p++)
    {
      for (lin = p; lin <= tamSeqMenor; lin += size - 1)
      {
        MPI_Recv(matrizEscores[lin], tamSeqMaior + 1, MPI_INT, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
  }
  else
  {
    // Envia as linhas calculadas para o processo 0
    for (lin = rank; lin <= tamSeqMenor; lin += size - 1)
    {
      MPI_Send(matrizEscores[lin], tamSeqMaior + 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
  }

  if (rank == 0)
  {

    // Definir o maior e o menor valor da matriz
    linPMaior = 1;
    colPMaior = 1;
    PMaior = matrizEscores[1][1];

    linUMaior = 1;
    colUMaior = 1;
    UMaior = matrizEscores[1][1];

    for (int lin = 1; lin <= tamSeqMenor; lin++)
    {
      for (int col = 1; col <= tamSeqMaior; col++)
      {
        if (PMaior < matrizEscores[lin][col])
        {
          linPMaior = lin;
          colPMaior = col;
          PMaior = matrizEscores[lin][col];
        }
        if (UMaior <= matrizEscores[lin][col])
        {
          linUMaior = lin;
          colUMaior = col;
          UMaior = matrizEscores[lin][col];
        }
      }
    }

    printf("\nMatriz de escores Gerada.");
    printf("\nPrimeiro Maior escore = %d na celula [%d,%d]", PMaior, linPMaior, colPMaior);
    printf("\nUltimo Maior escore = %d na celula [%d,%d]", UMaior, linUMaior, colUMaior);
  }
}
/* imprime a matriz de escores de acordo */
void mostraMatrizEscores()
{
  int i, lin, col;

  printf("\nMatriz de escores Atual:\n");

  printf("%4c%4c", ' ', ' ');
  for (i = 0; i <= tamSeqMaior; i++)
    printf("%4d", i);
  printf("\n");

  printf("%4c%4c%4c", ' ', ' ', '-');
  for (i = 0; i < tamSeqMaior; i++)
    printf("%4c", mapaBases[(seqMaior[i])]);
  printf("\n");

  printf("%4c%4c", '0', '-');
  for (col = 0; col <= tamSeqMaior; col++)
    printf("%4d", matrizEscores[0][col]);
  printf("\n");

  for (lin = 1; lin <= tamSeqMenor; lin++)
  {
    printf("%4d%4c", lin, mapaBases[(seqMenor[lin - 1])]);
    for (col = 0; col <= tamSeqMaior; col++)
    {
      printf("%4d", matrizEscores[lin][col]);
    }
    printf("\n");
  }
}

/* mostra os alinhamentos */
void mostraAlinhamentoGlobal(void)
{
  int i;

  printf("\nAlinhamento Obtido - Tamanho = %d:\n", tamAlinha);

  printf("%c", mapaBases[alinhaGMaior[0]]);
  for (i = 1; i < tamAlinha; i++)
    printf("%c", mapaBases[alinhaGMaior[i]]);
  printf("\n");

  printf("%c", mapaBases[alinhaGMenor[0]]);
  for (i = 1; i < tamAlinha; i++)
    printf("%c", mapaBases[alinhaGMenor[i]]);
  printf("\n");
}

void salvaMatrizEmArquivo(const char *nomeArquivo)
{
  FILE *arquivo = fopen(nomeArquivo, "w");
  if (arquivo == NULL)
  {
    perror("Erro ao abrir o arquivo");
    exit(EXIT_FAILURE);
  }

  fprintf(arquivo, "Matriz de escores Atual:\n");

  fprintf(arquivo, "%4c%4c", ' ', ' ');
  for (int i = 0; i <= tamSeqMaior; i++)
  {
    fprintf(arquivo, "%4d", i);
  }
  fprintf(arquivo, "\n");

  fprintf(arquivo, "%4c%4c%4c", ' ', ' ', '-');
  for (int i = 0; i < tamSeqMaior; i++)
  {
    fprintf(arquivo, "%4c", mapaBases[(seqMaior[i])]);
  }
  fprintf(arquivo, "\n");

  fprintf(arquivo, "%4c%4c", '0', '-');
  for (int col = 0; col <= tamSeqMaior; col++)
  {
    fprintf(arquivo, "%4d", matrizEscores[0][col]);
  }
  fprintf(arquivo, "\n");

  for (int lin = 1; lin <= tamSeqMenor; lin++)
  {
    fprintf(arquivo, "%4d%4c", lin, mapaBases[(seqMenor[lin - 1])]);
    for (int col = 0; col <= tamSeqMaior; col++)
    {
      fprintf(arquivo, "%4d", matrizEscores[lin][col]);
    }
    fprintf(arquivo, "\n");
  }

  fclose(arquivo);
  printf("Matriz de scores salva no arquivo '%s'\n", nomeArquivo);
}

/* gera o alinhamento global por meio do percurso de retorno na Matriz de escores,
   em duas formas possiveis, conforme o parametro tipo: 1) a partir da celula do
   primeiro maior escore [linPMaior,colPMaior] ou 2) a partir da celula de ultimo
   maior escore [linUMaior,colUMaior], ambos em direcao a celula inicial [0,0],
   uma celula por vez. O caminho de retorno deve ser feito seguindo o mesmo caminho
   inverso que gerou a celular final a partir da celula inicial. O alinhamento
   global eh composto por duas sequencias alinhaGMenor e alinhaGMaior.

   Note que outros alinhamentos globais podem ser encontrados se o percurso do
   traceback for ramificado em celulas que foram "escoreadas/pontuadas" por meio
   de uma estrategia de desempate, pelo fato de ter havido empate na pontuacao
   dos nohs vizinhos. Alem disso, alinhamentos parciais tambem podem ser obtidos
   com traceback iniciado a partir de qualquer c�lula */
void traceBack(int tipo)
{
  int tbLin, tbCol, peso, pos, aux, i;

  // O processo 0 deve ser o único a realizar o traceback
  if (tipo == 1)
  {
    printf("\nGeracao do Primeiro Maior Alinhamento Global:\n");
    tbLin = linPMaior;
    tbCol = colPMaior;
  }
  else
  {
    printf("\nGeracao do Ultimo Maior Alinhamento Global:\n");
    tbLin = linUMaior;
    tbCol = colUMaior;
  }

  pos = 0;
  do
  {
    if (tbLin > 0 && tbCol > 0)
    {
      // Verifica o escore do elemento [tbLin, tbCol]
      peso = matrizPesos[(seqMenor[tbLin - 1])][(seqMaior[tbCol - 1])];
      escoreDiag = matrizEscores[tbLin - 1][tbCol - 1] + peso;
      escoreLin = matrizEscores[tbLin][tbCol - 1] - penalGap;
      escoreCol = matrizEscores[tbLin - 1][tbCol] - penalGap;

      if ((escoreDiag > escoreLin) && (escoreDiag > escoreCol))
      {
        // Se houver um gap duplo
        if (seqMenor[tbLin - 1] != seqMaior[tbCol - 1])
        {
          printf("\nALERTA no TraceBack: Pos = %d Lin = %d e Col = %d\n", pos, tbLin, tbCol);

          alinhaGMenor[pos] = X;
          alinhaGMaior[pos] = seqMaior[tbCol - 1];
          tbCol--;
          pos++;

          alinhaGMenor[pos] = seqMenor[tbLin - 1];
          alinhaGMaior[pos] = X;
          tbLin--;
          pos++;
        }
        else
        {
          alinhaGMenor[pos] = seqMenor[tbLin - 1];
          tbLin--;
          alinhaGMaior[pos] = seqMaior[tbCol - 1];
          tbCol--;
          pos++;
        }
      }
      else if (escoreLin >= escoreCol)
      {
        alinhaGMenor[pos] = X;
        alinhaGMaior[pos] = seqMaior[tbCol - 1];
        tbCol--;
        pos++;
      }
      else
      {
        alinhaGMenor[pos] = seqMenor[tbLin - 1];
        alinhaGMaior[pos] = X;
        tbLin--;
        pos++;
      }
    }
    else
    {
      break; // Se tbLin ou tbCol forem <= 0, saímos do loop para evitar acessos inválidos
    }

  } while ((tbLin != 0) && (tbCol != 0));

  /* descarrega o restante de gaps da linha 0, se for o caso */
  while (tbLin > 0)
  {
    alinhaGMenor[pos] = seqMenor[tbLin - 1];
    alinhaGMaior[pos] = X;
    tbLin--;
    pos++;
  }

  /* descarrega o restante de gaps da coluna 0, se for o caso */
  while (tbCol > 0)
  {
    alinhaGMenor[pos] = X;
    alinhaGMaior[pos] = seqMaior[tbCol - 1];
    tbCol--;
    pos++;
  }

  tamAlinha = pos;

  /* Inverte o alinhamento para corrigir a ordem */
  for (i = 0; i < (tamAlinha / 2); i++)
  {
    aux = alinhaGMenor[i];
    alinhaGMenor[i] = alinhaGMenor[tamAlinha - i - 1];
    alinhaGMenor[tamAlinha - i - 1] = aux;

    aux = alinhaGMaior[i];
    alinhaGMaior[i] = alinhaGMaior[tamAlinha - i - 1];
    alinhaGMaior[tamAlinha - i - 1] = aux;
  }

  printf("\nAlinhamento Global Gerado.");
}

/* menu de opcoes fornecido para o usuario */
int menuOpcao(void)
{
  int op;
  char enter;

  do
  {
    printf("\nMenu de Opcao:");
    printf("\n<01> Ler Matriz de Pesos");
    printf("\n<02> Mostrar Matriz de Pesos");
    printf("\n<03> Ler Penalidade de Gap");
    printf("\n<04> Mostrar Penalidade");
    printf("\n<05> Definir Sequencias Genomicas");
    printf("\n<06> Mostrar Sequencias");
    printf("\n<07> Gerar Matriz de Escores");
    printf("\n<08> Mostrar Matriz de Escores");
    printf("\n<09> Gerar Alinhamento Global");
    printf("\n<10> Mostrar Alinhamento Global");
    printf("\n<11> Sair");
    printf("\nDigite a opcao => ");
    scanf("%d", &op);
    scanf("%c", &enter);
  } while ((op < 1) || (op > sair));

  return (op);
}

/* trata a opcao fornecida pelo usuario, executando o modulo pertinente */
void trataOpcao(int op, int rank, int size)
{
  int resp;
  char enter;
  char fileName[100];

  switch (op)
  {
  case 1:
    leMatrizPesos(rank);
    break;
  case 2:
    if (rank == 0)
      mostraMatrizPesos();
    break;
  case 3:
    penalGap = lePenalidade(rank);
    MPI_Bcast(&penalGap, 1, MPI_INT, 0, MPI_COMM_WORLD);
    break;
  case 4:
    printf("\nPenalidade = %d", penalGap);
    break;
  case 5:
    printf("\nDeseja Definicao: <1>MANUAL, <2>ALEATORIA? ou <3>ARQUIVO= ");
    scanf("%d", &resp);
    scanf("%c", &enter); /* remove o enter */

    resp_geracao = resp;
    MPI_Bcast(&resp_geracao, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (resp == 1)
    {
      leSequencias(rank);
    }
    if (resp == 2)
    {
      leTamMaior(rank);
      leTamMenor(rank);
      grauMuta = leGrauMutacao(rank);
      geraSequencias(rank);
    }
    if (resp == 3)
    {
      printf("Digite o nome do arquivo das sequencias: ");
      scanf("%s", fileName);
      leSequenciasDeArquivo(fileName, rank);
    }
    break;
  case 6:
    if (rank == 0)
      mostraSequencias();
    break;
  case 7:
    // Processo 0 pode salvar a matriz e fazer o traceback
    geraMatrizEscores(rank, size); // Chamada da versão paralela
    if (rank == 0)
    {
      printf("\nGerando matriz de escores em paralelo com MPI...\n");
      salvaMatrizEmArquivo("matriz_escores.txt");
    }
    break;
  case 8:
    if (rank == 0)
      mostraMatrizEscores();
    break;
  case 9:
    if (rank == 0)
    {
      printf("\nDeseja: <1> Primeiro Maior ou <2> Ultimo Maior? = ");
      scanf("%d", &resp);
      scanf("%c", &enter); /* remove o enter */
      traceBack(resp);
    }
    break;
  case 10:
    if (rank == 0)
      mostraAlinhamentoGlobal();
    break;
  }
}
/* programa principal */
void main(int argc, char *argv[])
{
  int opcao;
  int rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  srand(time(NULL));

  if (rank == 0)
  {
    do
    {
      printf("\n\nPrograma Needleman-Wunsch Paralelo\n");
      opcao = menuOpcao();

      // Broadcast da opção para todos os processos
      MPI_Bcast(&opcao, 1, MPI_INT, 0, MPI_COMM_WORLD);

      // Processo 0 gerencia as opções, mas não faz o cálculo
      trataOpcao(opcao, rank, size);

    } while (opcao != sair);
  }
  else
  {
    while (1)
    {
      // Recebe o broadcast da opção
      MPI_Bcast(&opcao, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if (opcao == 5)
      {
        MPI_Bcast(&resp_geracao, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (resp_geracao == 1)
        {
          MPI_Bcast(&tamSeqMaior, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(&tamSeqMenor, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(seqMaior, tamSeqMaior, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(seqMenor, tamSeqMenor, MPI_INT, 0, MPI_COMM_WORLD);

          printf("Processo %d", rank);
          printf("\n");

          for (int i = 0; i < tamSeqMaior; i++)
          {
            printf("%d ", seqMaior[i]);
          }

          printf("\n");

          for (int i = 0; i < tamSeqMenor; i++)
          {
            printf("%d ", seqMenor[i]);
          }

          printf("\n");
          continue;
        }
        if (resp_geracao == 2)
        {
          MPI_Bcast(&prob, 1, MPI_INT, 0, MPI_COMM_WORLD);

          MPI_Bcast(&tamSeqMaior, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(&tamSeqMenor, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(seqMaior, tamSeqMaior, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(seqMenor, tamSeqMenor, MPI_INT, 0, MPI_COMM_WORLD);
          printf("Processo %d", rank);
          printf("\n");

          for (int i = 0; i < tamSeqMaior; i++)
          {
            printf("%d ", seqMaior[i]);
          }

          printf("\n");

          for (int i = 0; i < tamSeqMenor; i++)
          {
            printf("%d ", seqMenor[i]);
          }

          printf("\n");
          continue;
        }
        if (resp_geracao == 3)
        {
          MPI_Bcast(&tamSeqMaior, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(&tamSeqMenor, 1, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(seqMaior, tamSeqMaior, MPI_INT, 0, MPI_COMM_WORLD);
          MPI_Bcast(seqMenor, tamSeqMenor, MPI_INT, 0, MPI_COMM_WORLD);

          printf("Processo %d", rank);
          printf("\n");

          for (int i = 0; i < tamSeqMaior; i++)
          {
            printf("%d ", seqMaior[i]);
          }

          printf("\n");

          for (int i = 0; i < tamSeqMenor; i++)
          {
            printf("%d ", seqMenor[i]);
          }

          printf("\n");
          continue;
        }
      }

      if (opcao != 9 && opcao != 10)
        trataOpcao(opcao, rank, size);
    }
  }

  MPI_Finalize();
}