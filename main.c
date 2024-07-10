/*

Programa de Alinhamento de Genoma Vers�o 28/03/2022.

Programa "sequencial" baseado no metodo Needleman-Wunsch, desenvolvido
por Saul B. Needleman e Christian D. Wunsch, publicado em 1970, conforme segue:

Needleman, Saul B. & Wunsch, Christian D. (1970).
"A general method applicable to the search for similarities in the
amino acid sequence of two proteins". Journal of Molecular Biology.
48 (3):443-53.

Este programa NAO pode ser usado, da forma original ou modificada, completa
ou parcial, para finalidade comercial ou lucrativa, sem a previa autorizacao
dos seus desenvolvedores, os quais detem os direitos autorais.

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
Data da modificacao: 28/03/2022
eMail: ralgonca@gmail.com
Whatsapp: 55(44)99159-2310
Modificacoes realizadas: desenvolvimento do codigo original, composto pelos modulos
leTamMaior(), leTamMenor(), lePenalidade(), menuOpcao(void), trataOpcao(), geraSequencias(),
leMatrizPesos(), mostraMatrizPesos(), leGrauMutacao(), geraMatrizScores(), mostraMatrizScores(),
leSequencias(), geraSequencias(), mostraSequencias(), traceBack(), mostraAlinhamentoGlobal()
e main().



======================================
BREVE DESCRI��O:

======================================

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define A 0
#define T 1
#define G 2
#define C 3
#define sair 11

#define maxSeq 1000

/* baseMapa mapeia indices em caracteres que representam as bases,
   sendo 0='A', 1='T', 2='G', 3='C' e 4='-' representando gap */
char baseMapa[5]={'A','T','G','C','-'};

/* seqMaior e seqMenor representam as duas sequencias de bases de
   entrada, a serem comparadas, inicializadas conforme segue. Elas
   conterao os indices aos inves dos proprios caracteres. seqMenor
   deve ser menor ou igual a seqMaior. */
char seqMaior[maxSeq]={A,A,C,T,T,A},
     seqMenor[maxSeq]={A,C,T,T,G,A},
     alinhaMaior[maxSeq],
     alinhaMenor[maxSeq];

/* matrizScores representa a matriz de scores que sera preenchida
   pelo metodo. A matriz, ao final de seu preenchimento, permitira
   obter o melhor alinhamento global entre as sequencias seqMaior e
   seqMenor, por meio de uma operacao denominada TraceBack. Uma linha
   e uma coluna extras sao adicionadas na matriz para correlacionar
   as bases com gaps. Trata-se da linha 0 e coluna 0. A matriz de scores
   tera tamSeqMenor+1 linhas e tamSeqMaior+1 colunas. Considera-se a
   primeira dimensao da matriz como linhas e a segunda como colunas.*/
int matrizScores[maxSeq+1][maxSeq+1];

int tamSeqMaior=6,  /* tamanho da sequencia maior, inicializado como 6 */
    tamSeqMenor=6,  /* tamanho da sequencia menor, inicializado como 6 */
    tamAlinha,      /* tamanho do alinhamento global obtido */
    penalGap=0,     /* penalidade de gap, a ser descontada no score acumulado */
    grauMuta=0,     /* porcentagem maxima de mutacao na geracao aleatoria
                       da sequencia menor, a qual eh copiada da maior e
                       sofre algumas trocas de bases */
    diagScore,      /* score da diagonal anterior da matriz de scores */
    linScore,       /* score da linha anterior da matriz de scores */
    colScore;       /* score da coluna anterior da matriz de scores */

/*  matrizPesos contem os pesos do pareamento de bases. Estruturada e
    inicializada conforme segue, onde cada linha ou coluna se refere a
    uma das bases A, T, G ou C. Considera-se a primeira dimensao da
    matriz como linhas e a segunda como colunas.

       A T G C
       0 1 2 3
   A 0 1 0 0 0
   T 1 0 1 0 0
   G 2 0 0 1 0
   C 3 0 0 0 1
*/
int matrizPesos[4][4]={1z,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};


/* leitura do tamanho da sequencia maior */
int leTamMaior(void)
{
  printf("\nLeitura do Tamanho da Sequencia Maior:");
  do
  { printf("\nDigite 0 < valor < %d = ", maxSeq);
    scanf("%d", &tamSeqMaior);
  } while ((tamSeqMaior<1)||(tamSeqMaior>maxSeq));
}

/* leitura do tamanho da sequencia menor */
int leTamMenor(void)
{
  printf("\nLeitura do Tamanho da Sequencia Menor:");
  do
  {  printf("\nDigite 0 < valor <= %d = ", tamSeqMaior);
     scanf("%d", &tamSeqMenor);
  } while ((tamSeqMenor<1)||(tamSeqMenor>tamSeqMaior));
}

/* leitura do valor da penalidade de gap */
int lePenalidade(void)
{   int penal;

  printf("\nLeitura da Penalidade de Gap:");
  do
  {
    printf("\nDigite valor >= 0 = ");
    scanf("%d", &penal);
  } while (penal<0);

  return penal;
}

/* leitura da matriz de pesos */
void leMatrizPesos()
{ int i,j;

  printf("\nLeitura da Matriz de Pesos:\n");
  for (i=0; i<4; i++)
  {
    for (j=0; j<4; j++)
    {
      printf("Digite valor %c x %c = ", baseMapa[i], baseMapa[j]);
      scanf("%d",&(matrizPesos[i][j]));
    }
    printf("\n");
  }
}

/* mostra da matriz de pesos */
void mostraMatrizPesos(void)
{ int i,j;

  printf("\nMatriz de Pesos Atual:");
  printf("\n%4c%4c%4c%4c%4c\n",' ','A','T','G','C');
  for (i=0; i<4; i++)
  {
    printf("%4c",baseMapa[i]);
    for (j=0; j<4; j++)
      printf("%4d",matrizPesos[i][j]);
    printf("\n");
  }
}


/* leitura da porcentagem maxima (grau) de mutacao aleatoria. Essa
   porcentagem eh usada na geracao aleatoria da seqMenor. A seqMenor,
   apohs extraida da seqMaior, sobre algumas trocas de bases, para
   causar alguma diferenciacao aleatoria. A quantidade de trocas
   realizadas eh no maximo a porcentagem aqui informada. */

int leGrauMutacao(void)
{ int prob;

  printf("\nLeitura da Porcentagem Maxima de Mutacao Aleatoria:\n");
  do
  { printf("\nDigite 0 <= valor <= 100 = ");
    scanf("%d", &prob);
  } while ((prob<0)||(prob>100));

  return prob;
}

/* leitura manual das sequencias de entrada seqMaior e seqMenor */
void leSequencias(void)
{ int i, erro;

  printf("\nLeitura das Sequencias:\n");

  /* lendo a sequencia maior */
  do
  {
    printf("\nPara a Sequencia Maior,");
    printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
    do
    { printf("\n> ");
      fgets(seqMaior,maxSeq,stdin);
      tamSeqMaior=strlen(seqMaior)-1; /* remove o enter */
    } while (tamSeqMaior<1);
    printf("\ntamSeqMaior = %d\n",tamSeqMaior);
    i=0;
    erro=0;
    do
    {
      switch (seqMaior[i])
      {
        case 'A': seqMaior[i]=(char)A;
                  break;
        case 'T': seqMaior[i]=(char)T;
                  break;
        case 'G': seqMaior[i]=(char)G;
                  break;
        case 'C': seqMaior[i]=(char)C;
                  break;
        default: erro=1;  /* nao eh permitido qquer outro caractere */
      }
      i++;
    } while ((erro==0)&&(i<tamSeqMaior));
  }while (erro==1);

  /* lendo a sequencia maior */
  do
  {
    printf("\nPara a Sequencia Menor, ");
    printf("\nDigite apenas caracteres 'A', 'T', 'G' e 'C'");
    do
    { printf("\n> ");
      fgets(seqMenor,maxSeq,stdin);
      tamSeqMenor=strlen(seqMenor)-1; /* remove o enter */
    } while ((tamSeqMenor<1)||(tamSeqMenor>tamSeqMaior));
    printf("\ntamSeqMenor = %d\n",tamSeqMenor);

    i=0;
    erro=0;
    do
    {
      switch (seqMenor[i])
      {
        case 'A': seqMenor[i]=(char)A;
                  break;
        case 'T': seqMenor[i]=(char)T;
                  break;
        case 'G': seqMenor[i]=(char)G;
                  break;
        case 'C': seqMenor[i]=(char)C;
                  break;
        default: erro=1;
      }
      i++;
    } while ((erro==0)&&(i<tamSeqMenor));
  }while (erro==1);
}


/* geracao das sequencias aleatorias, conforme tamanho. Gera-se numeros
   aleatorios de 0 a 3 representando as bases 'A', 'T', 'G' e 'C'.
   Gera-se primeiramente a maior sequencia e desta extrai a menor
   sequencia. A menor sequencia ainda sofre algumas trocas de bases ou
   mutacoes, de acordo com o grau de mutacao informado. A ideia eh
   gerar sequencias parecidas, mas com certo grau de diferenca. */

void geraSequencias(void)
{   int i, dif, probAux, ind, nTrocas;
    char base;

    srand(time(NULL));

    printf("\nGeracao Aleatoria das Sequencias:\n");

    /* gerando a sequencia maior */
    for (i=0; i<tamSeqMaior; i++)
      {
        base=(char)(rand()%4); /* produz valores de 0 a 3 e os converte
                                  em char */
        seqMaior[i]= base;
      }

    dif=tamSeqMaior-tamSeqMenor; /* diferenca entre os tamanhos das sequencias */

    ind=0;
    if (dif>0)
      ind=rand()%dif; /* produz um indice aleatorio na sequencia maior
                       para a partir dela copiar a sequencia menor */

    /* gerando a sequencia menor a partir da maior. Copia trecho da
       sequencia maior, a partir de um indice aleatorio que nao
       ultrapasse os limites do vetor maior */
    for (i=0; i<tamSeqMenor; i++)
        seqMenor[i]=seqMaior[ind+i];

    /* causa mutacoes aleatorias na sequencia menor para gerar "gaps",
       sobre cada base, de acordo com o grau (porcentagem) informado.
       A mutacao causa a troca da base original por outra base aleatoria
       necessariamente diferente. Gera-se uma probabilidade aleatoria
       ateh 100 e se ela estiver dentro do grau informado, a mutacao
       ocorre na base, caso contrario, mantem a base copiada. */

    i=0;
    nTrocas=0;
    while ((i<tamSeqMenor)&&(nTrocas<grauMuta))
    {
      probAux=rand()%100;

      if (probAux<grauMuta)
      {
        seqMenor[i]=(seqMenor[i]+(rand()%3)+1)%4;
        nTrocas++;
      }
    }

    printf("\nSequencias Geradas, Dif = %d Ind = %d\n",dif, ind);
}

/* mostra das sequencias seqMaior e seqMenor */
void mostraSequencias(void)
{   int i;

  printf("\nSequencias Atuais:\n");
  printf("\nSequencia Maior, Tam = %d\n", tamSeqMaior);
  printf("%c",baseMapa[(int)seqMaior[0]]);
  for (i=1; i<tamSeqMaior; i++)
    printf("%c",baseMapa[(int)seqMaior[i]]);
  printf("\n");

  printf("\nSequencia Menor, Tam = %d\n", tamSeqMenor);
  printf("%c",baseMapa[(int)seqMenor[0]]);
  for (i=1; i<tamSeqMenor; i++)
    printf("%c",baseMapa[(int)seqMenor[i]]);
  printf("\n");
}

/* geraMatrizScores gera a matriz de scores. A matriz de scores tera
   tamSeqMenor+1 linhas e tamSeqMaior+1 colunas. A linha 0 e a coluna
   0 s�o adicionadas para representar gaps e conter penalidades. As
   demais linhas e colunas s�o associadas as bases da seqMenor e da
   SeqMaior, respectivamente. */

void geraMatrizScores(void) /* TODO - fazer paralelismo aqui */
{ int lin, col, peso, linMaior, colMaior, maior;

  printf("\nGeracao da Matriz de Scores:\n");
  /*  A matriz ser� gerada considerando que ela representa o cruzamento
      da seqMenor[] associada as linhas e a seqMaior[] associada as
      colunas. */

  /* inicializando a linha de penalidades/gaps */
  for (col=0; col<=tamSeqMaior; col++)
    matrizScores[0][col]=-1*(col*penalGap);

  /* inicializando a coluna de penalidades/gaps */
  for (lin=0; lin<=tamSeqMenor; lin++)
    matrizScores[lin][0]=-1*(lin*penalGap);

  /* calculando os demais scores, percorrendo todas as posicoes
     da matriz, linha por linha, coluna por coluna, aplicando
     a seguinte f�rmula:
                             / f(lin-1,col-1)+matrizPesos[lin,col]
     f(lin,col) = m�ximo de {  f(lin,col-1)-penalGap
                             \ f(lin-1,col)-penalGap
  */

  for (lin=1; lin<=tamSeqMenor; lin++)
  {
    for (col=1; col<=tamSeqMaior; col++)
    {
      peso=matrizPesos[(int)(seqMenor[lin-1])][(int)(seqMaior[col-1])];
      diagScore = matrizScores[lin-1][col-1]+peso;
      linScore = matrizScores[lin][col-1]-penalGap;
      colScore = matrizScores[lin-1][col]-penalGap;

      if ((diagScore>=linScore)&&(diagScore>=colScore))
        matrizScores[lin][col]=diagScore;
      else if (linScore>colScore)
              matrizScores[lin][col]=linScore;
           else matrizScores[lin][col]=colScore;
    }
  }

  /* localiza o ultimo maior score, e sua posicao */

  linMaior=1;
  colMaior=1;
  maior=matrizScores[1][1];
  for (lin=1; lin<=tamSeqMenor; lin++)
  {
    for (col=1; col<=tamSeqMaior; col++)
    {
      if (maior<=matrizScores[lin][col])
      {
        linMaior=lin;
        colMaior=col;
        maior=matrizScores[lin][col];
      }
    }
  }
  printf("\nMatriz de Scores Gerada.");
  printf("\nUltimo Maior Score = %d na celula [%d,%d]", maior, linMaior, colMaior);
}


/* imprime a matriz de scores de acordo */
void mostraMatrizScores()
{ int i, lin, col;

  printf("\nMatriz de Scores Atual:\n");

  printf("%4c%4c",' ',' ');
  for (i=0; i<=tamSeqMaior; i++)
    printf("%4d",i);
  printf("\n");

  printf("%4c%4c%4c",' ',' ','-');
  for (i=0; i<tamSeqMaior; i++)
    printf("%4c",baseMapa[(int)(seqMaior[i])]);
  printf("\n");

  printf("%4c%4c",'0','-');
  for (col=0; col<=tamSeqMaior; col++)
    printf("%4d",matrizScores[0][col]);
  printf("\n");

  for (lin=1;lin<=tamSeqMenor;lin++)
  {
    printf("%4d%4c",lin,baseMapa[(int)(seqMenor[lin-1])]);
    for (col=0;col<=tamSeqMaior;col++)
    {
      printf("%4d",matrizScores[lin][col]);
    }
    printf("\n");
  }
}


/* mostra os alinhamentos */
void mostraAlinhamentoGlobal(void)
{   int i;

  printf("\nAlinhamentos Atuais - Tamanho = %d:\n", tamAlinha);

  printf("%c",baseMapa[(int)alinhaMaior[0]]);
  for (i=1; i<tamAlinha; i++)
    printf("%c",baseMapa[(int)alinhaMaior[i]]);
  printf("\n");

  printf("%c",baseMapa[(int)alinhaMenor[0]]);
  for (i=1; i<tamAlinha; i++)
    printf("%c",baseMapa[(int)alinhaMenor[i]]);
  printf("\n");
}

/* gera o alinhamento global por meio do retorno na Matriz de Scores
   da celula final (tamSeqMenor,tamSeqMaior) ateh a celula inicial (0,0).
   A partir da celula final, retorna-se para a celula de onde o score
   atual foi derivado. Repete-se esse processo a partir desta celula
   retornada e assim sucessivamente, ateh alcancar a celula inicial.
   O alinhamento global � composto por duas sequencias alinhaMenor e
   alinhaMaior. */
void traceBack() /* TODO - fazer paralelismo aqui */
{ int tbLin, tbCol, peso, pos, posAux, aux, i;

  printf("\nGeracao do Alinhamento Global:\n");
  tbLin=tamSeqMenor;
  tbCol=tamSeqMaior;
  pos=0;
  do
  {
    peso=matrizPesos[(int)(seqMenor[tbLin-1])][(int)(seqMaior[tbCol-1])];
    diagScore = matrizScores[tbLin-1][tbCol-1]+peso;
    linScore = matrizScores[tbLin][tbCol-1]-penalGap;
    colScore = matrizScores[tbLin-1][tbCol]-penalGap;

      if ((diagScore>=linScore)&&(diagScore>=colScore))
      {
        alinhaMenor[pos]=seqMenor[tbLin-1];
        alinhaMaior[pos]=seqMaior[tbCol-1];
        tbLin--;
        tbCol--;
        pos++;
      }
      else if (linScore>colScore)
           {
              alinhaMenor[pos]=(char)4;
              alinhaMaior[pos]=seqMaior[tbCol-1];
              tbCol--;
              pos++;
           }
           else
           {
              alinhaMenor[pos]=seqMenor[tbLin-1];
              alinhaMaior[pos]=(char)4;
              tbLin--;
              pos++;
           }

  } while ((tbLin!=0)&&(tbCol!=0));

  /* descarrega o restante de gaps da linha 0, se for o caso */
  while (tbLin>0)
  {
    alinhaMenor[pos]=seqMenor[tbLin-1];
    alinhaMaior[pos]=(char)4;
    tbLin--;
    pos++;
  }

  /* descarrega o restante de gaps da coluna 0, se for o caso */
  while (tbCol>0)
  {
    alinhaMenor[pos]=(char)4;
    alinhaMaior[pos]=seqMaior[tbCol-1];
    tbCol--;
    pos++;
  }

  tamAlinha=pos;

  /* o alinhamento foi feito de tras para frente e deve ser
     invertido, conforme segue */
  for (i=0;i<(tamAlinha/2);i++)
  {
    aux=alinhaMenor[i];
    alinhaMenor[i]=alinhaMenor[tamAlinha-i-1];
    alinhaMenor[tamAlinha-i-1]=aux;

    aux=alinhaMaior[i];
    alinhaMaior[i]=alinhaMaior[tamAlinha-i-1];
    alinhaMaior[tamAlinha-i-1]=aux;
  }

  printf("\nAlinhamento Global Gerado.");
}

/* menu de opcoes fornecido para o usuario */
int menuOpcao(void)
{ int op;
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
    printf("\n<07> Gerar Matriz de Scores");
    printf("\n<08> Mostrar Matriz de Scores");
    printf("\n<09> Gerar Alinhamento Global");
    printf("\n<10> Mostrar Alinhamento Global");
    printf("\n<11> Sair");
    printf("\nDigite a opcao => ");
    scanf("%d",&op);
    scanf("%c",&enter);
  } while ((op<1)||(op>sair));

  return (op);
}

/* trata a opcao fornecida pelo usuario, executando o modulo pertinente */
void trataOpcao(int op)
{ int resp;
  char enter;

  switch (op)
  {
    case 1: leMatrizPesos();
            break;
    case 2: mostraMatrizPesos();
            break;
    case 3: penalGap=lePenalidade();
            break;
    case 4: printf("\nPenalidade = %d",penalGap);
            break;
    case 5: printf("\nDeseja Definicao: <1>MANUAL ou <2>ALEATORIA? = "); /* TODO - leitura do arquivo deve ser inserida aqui */
            scanf("%d",&resp);
            scanf("%c",&enter); /* remove o enter */
            if (resp==1)
            {
              leSequencias();
            }
            else
            { leTamMaior();
              leTamMenor();
              grauMuta=leGrauMutacao();
              geraSequencias();
            }
            break;
    case 6: mostraSequencias();
            break;
    case 7: geraMatrizScores();
            break;
    case 8: mostraMatrizScores();
            break;
    case 9: traceBack();
            break;
    case 10: mostraAlinhamentoGlobal();
            break;
  }
}

/* programa principal */
void main(void)
{ int opcao;

  do
  {
    printf("\n\nPrograma Needleman-Wunsch Sequencial\n");
    opcao=menuOpcao();
    trataOpcao(opcao);

  } while (opcao!=sair);

}
