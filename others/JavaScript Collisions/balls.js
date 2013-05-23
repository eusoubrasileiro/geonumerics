var STEP = 1;
var NUMBERBALLS = 100;
var SIZE = 500.0; /* Posicao máximo em x e y */
var RADIUS = 1.0; /* Raio das bolinhas */
var MASS = 1.0;
var ERROR = 0.5; /* Erro de ultrapassagem dos limites de ambiente */

/* collection of balls */
var BallSet = [];

function randsign()
{
	return (Math.random() > 0.5) ? -1.0 : 1.0;
}

function randomcolor() {
    var letters = '0123456789ABCDEF'.split('');
    var color = '#';
    for (var i = 0; i < 6; i++ ) {
        color += letters[Math.round(Math.random() * 15)];
    }
    return color;
}

function Ball(x, y, vi, vj, mass, radius, outx, outy, ignore)
{
	this.x = x;
	this.y = y;
	this.vi = vi;
	this.vj = vj;
	this.mass = mass;
	this.radius = radius;
	/* indicativo de dentro (0) ou fora do ambiente (1) */
	this.outx = outx;
	this.outy = outy;
	this.ignore = ignore;
	this.color = randomcolor();

	this.Move = function()
	{ /* move in the versor direction one STEP*/
		this.x += this.vi*STEP;
		this.y += this.vj*STEP;
	}
}


function Start(size, numberballs, step, radius, mass)
{    

    SIZE = size;
    NUMBERBALLS = numberballs;
    STEP = step;
    MASS = mass;
    RADIUS = radius;

    for(var i=0; i<NUMBERBALLS; i++){

		/* posicoes entre 0 e SIZE para cada dimensao pois o ambiente its a cube (x, y) */
		/* after velocities maximum 1, mass and radius all equal */
    	var ball = new Ball(
    		(SIZE-2*RADIUS)*Math.random()+2*RADIUS, 
    		(SIZE-2*RADIUS)*Math.random()+2*RADIUS, 
    		randsign()*Math.random(),
    		randsign()*Math.random(),
    		MASS,
    		RADIUS,
    		false, false, false /* inside and not ignored */
    		);

    	BallSet.push(ball);
    }
}

/*
Loop infinito de
    (1) mover bola[i]
    (2) verificar colisao com ambiente
        (2.1) Caso haja -> inverte vetor velocidade
    (3) verficar colisao com todas outras bolas
        (3.1) Caso haja -> calcula novas velocidades por colisao elastica
*/
function Simulate(){
    
    for(var i=0; i<NUMBERBALLS; i++) /* limpa os ignores pra a nova rodada */
        BallSet[i].ignore = false;

    for(var i=0; i<NUMBERBALLS; i++){
    	BallSet[i].Move()        
        Collisions(i);
        CollisionBorders(i); /* da particula i */
    }
}

/*
compara a posicao do centro da bola i com o todas as demais j 
caso estejam mais proximao que a soma de seus raios COLISAO!!!
*/
function Collisions(i){
    
	for(var j=0; j < NUMBERBALLS; j++){	    
        if(!BallSet[j].ignore){ /* verfica se a colisao ja foi analisada com esta bola*/
            var inx = BallSet[i].x - BallSet[j].x;
            var iny = BallSet[i].y - BallSet[j].y;	            
            var dradius = BallSet[i].radius + BallSet[j].radius;
			/* verifica se as duas bolas estao mais proximas que a soma dos dois raios?? */
            if(inx*inx + iny*iny <= dradius*dradius)
                Collision(i, j)
            
        }
    }

    BallSet[i].ignore = true; /* apos calcular a colisao com todas ignorar*/
}

/*
COLISAO COMPLETAMENTE ELASTICA:::
conservacao de momento e conservacao de energia
em cada uma das direcoes huahuauha
pag. 172 - (9.4.10)
H. Moyses 1 - Física Básica - Mecânica
*/

function Collision(i, j){
    /* mr  mass reduzida razao entre mu e mi */
    /* p1f, p2f, p1i, p2i momentos finais de i (1) e j (2) */

    var mr = BallSet[j].mass/BallSet[i].mass;

    var p1i = BallSet[i].mass*BallSet[i].vi; /* alfa */
    var p2i = BallSet[j].mass*BallSet[j].vi; /* beta */
    /* conservacao de tudo em i/x para a bola[i]*/
    var p1f = (-p1i*(1.0-mr)/(1.0+mr))+(p2i*mr*2.0/(1.0+mr));/* phi */
    var p2f = (p1i*2.0/(1.0+mr))+(p2i*(1.0-mr)/(1.0+mr)); /* teta */
    BallSet[i].vi = p1f/BallSet[i].mass;
    BallSet[j].vi = p2f/BallSet[j].mass;

    p1i = BallSet[i].mass*BallSet[i].vj; /* alfa */
    p2i = BallSet[j].mass*BallSet[j].vj; /* beta */
    /* conservacao de tudo em j/y para a bola[i]*/
    p1f = (-p1i*(1.0-mr)/(1.0+mr))+(p2i*mr*2.0/(1.0+mr));/* phi */
    p2f = (p1i*2.0/(1.0+mr))+(p2i*(1.0-mr)/(1.0+mr)); /* teta */
    BallSet[i].vj = p1f/BallSet[i].mass;
    BallSet[j].vj = p2f/BallSet[j].mass;
}

/*
*) Se a bola esta dentro do ambiente
     *) Verifica se a bola saiu ou esta no limiar de sair
          *) Se sim inverte o vetor velocidade na direcao que ela saiu
     (* Se nao saiu nao faz nada
(* Se nao esta dentro nao faz nada - espera ela voltar
*/

function CollisionBorders(i){
    if(!BallSet[i].outx) /* se estava dentro pode ser que agora nao mais */
        if( BallSet[i].x+BallSet[i].radius>=SIZE || BallSet[i].x-BallSet[i].radius<=0 ){ /* verifica nesta dimensao */
            BallSet[i].vi *= -1; /* se saiu nesta dimensao inverte a velocidade */
            BallSet[i].outx = true; /* agora esta fora do ambiente nesta dimensao */
        }
    if(!BallSet[i].outy)
        if( BallSet[i].y+BallSet[i].radius >=SIZE|| BallSet[i].y-BallSet[i].radius<=0 ){ /* idem ... */
            BallSet[i].vj *= -1;
            BallSet[i].outy = true;
        }
    /* para verificar se ela jah voltou em alguma
    das tres dimensoes e agora estah dentro */
    if(BallSet[i].x<SIZE&&BallSet[i].x>0)
        BallSet[i].outx = false; /* se estah dentro entaum ja voltou logo out = 0 */
    if(BallSet[i].y<SIZE&&BallSet[i].y>0)
        BallSet[i].outy = false;   
}