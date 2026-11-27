#!/usr/bin/env python3
# deps: none (stdlib only; needs gcc or cc on PATH). Run: python3 ours/s20581.py [--hours 5.5] [--jobs N]
# L-meta-20581 search: SA over OR_k(P_k AND Q_k) cross-pair forms (P,Q mixed OR-sets of x/y literals) with
# Paar-shared OR chains, then CGP gate-count polishing of the best circuit. Writes ours/results/20581_best.json.
import os, sys, json, time, random, subprocess, tempfile, argparse, shutil
HERE = os.path.dirname(os.path.abspath(__file__))
ROWS = [0xf2bc7573, 0x59092a9d, 0xc443c8c2, 0x1d58483a, 0x1b7982bb, 0xa8882290, 0x74d0fd25, 0x00041402, 0x6879a8a5, 0xd65a6f1a, 0x7455fba3, 0x24508c1a, 0x3050a010, 0x68658ad5, 0x10a00735, 0x5e55690e, 0x6ab804a4, 0x00b01458, 0xf447af40, 0x00b07351, 0xa88c2290, 0x40d9d7bc, 0x40241486, 0x54058082, 0x6891afb9, 0x50ddf79e, 0x553080be, 0x28c180f8, 0x6ac8ae9e, 0xa6860e80, 0x20508c18, 0x10f01050]
PSA_C = r'''// SA over "OR of (P AND Q)" forms, P,Q mixed OR-sets of x/y literals; cost via Paar SLP greedy.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
typedef uint64_t u64; typedef uint32_t u32;
static u32 ROWS[32]={0xf2bc7573, 0x59092a9d, 0xc443c8c2, 0x1d58483a, 0x1b7982bb, 0xa8882290, 0x74d0fd25, 0x00041402, 0x6879a8a5, 0xd65a6f1a, 0x7455fba3, 0x24508c1a, 0x3050a010, 0x68658ad5, 0x10a00735, 0x5e55690e, 0x6ab804a4, 0x00b01458, 0xf447af40, 0x00b07351, 0xa88c2290, 0x40d9d7bc, 0x40241486, 0x54058082, 0x6891afb9, 0x50ddf79e, 0x553080be, 0x28c180f8, 0x6ac8ae9e, 0xa6860e80, 0x20508c18, 0x10f01050};
static u32 COLS[32];
#define MK 48
static u64 P[MK],Q[MK]; // bits 0..31 rows (x), 32..63 cols (y)
static int cnt[32][32]; static int unc;
static u64 rs=88172645463325252ULL;
static inline u64 rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static inline double rndf(void){ return (rnd()>>11)*(1.0/9007199254740992.0); }
static inline u32 rows_(u64 s){ return (u32)s; } static inline u32 cols_(u64 s){ return (u32)(s>>32); }
static int active_(int k){ return P[k] && Q[k]; }
// rect cells add
static void rect(u32 r,u32 c,int d){ if(!r||!c) return;
  for(int u=0;u<32;u++) if(r>>u&1) for(int v=0;v<32;v++) if(c>>v&1){ if(d>0){ if(cnt[u][v]++==0) unc--; } else { if(--cnt[u][v]==0) unc++; } } }
static void term(int k,int d){ if(!active_(k)) return; rect(rows_(P[k]),cols_(Q[k]),d); rect(rows_(Q[k]),cols_(P[k]),d); }
static int valid(u64 p,u64 q){
  if(rows_(p)&rows_(q)) return 0; if(cols_(p)&cols_(q)) return 0;
  u32 r=rows_(p), c=cols_(q); for(int u=0;u<32;u++) if(r>>u&1) if((ROWS[u]&c)!=c) return 0;
  r=rows_(q); c=cols_(p); for(int u=0;u<32;u++) if(r>>u&1) if((ROWS[u]&c)!=c) return 0;
  return 1;
}
// ---- Paar SLP cost ----
#define MW 1024
static unsigned short pc2[MW][MW];
static int paar(u64 *sets,int ns,int randomize,int (*rec)[3]){
  // dedupe & drop singletons
  static int el[2*MK][64]; static int len[2*MK]; int m=0;
  for(int i=0;i<ns;i++){ u64 s=sets[i]; if(__builtin_popcountll(s)<2) continue; int dup=0; for(int j=0;j<i;j++) if(sets[j]==s){dup=1;break;} if(dup) continue;
    len[m]=0; for(int b=0;b<64;b++) if(s>>b&1) el[m][len[m]++]=b; m++; }
  int nw=64, cost=0;
  static int tp[65536][2];
  while(1){
    int ntp=0; int best=0, ba=-1, bb=-1, nties=0;
    for(int i=0;i<m;i++){ if(len[i]<2) continue; for(int x=0;x<len[i];x++) for(int y=x+1;y<len[i];y++){ int a=el[i][x], b=el[i][y]; if(a>b){int t=a;a=b;b=t;}
        if(pc2[a][b]==0){ if(ntp<65536){ tp[ntp][0]=a; tp[ntp][1]=b; ntp++; } }
        int c=++pc2[a][b];
        if(c>best){ best=c; ba=a; bb=b; nties=1; } else if(c==best && randomize){ nties++; if(rnd()%nties==0){ba=a;bb=b;} } } }
    for(int t=0;t<ntp;t++) pc2[tp[t][0]][tp[t][1]]=0;
    if(best==0) break;
    if(nw>=MW) return 9999;
    int nwid=nw++; cost++;
    if(rec){ rec[cost-1][0]=ba; rec[cost-1][1]=bb; rec[cost-1][2]=nwid; }
    for(int i=0;i<m;i++){ int ha=-1,hb=-1; for(int x=0;x<len[i];x++){ if(el[i][x]==ba) ha=x; if(el[i][x]==bb) hb=x; }
      if(ha>=0&&hb>=0){ el[i][ha]=nwid; el[i][hb]=el[i][--len[i]]; } }
  }
  return cost;
}
static double LAM=3.0;
static int nact(void){ int n=0; for(int k=0;k<MK;k++) if(active_(k)) n++; return n; }
static int slpcost(int randomize){ u64 s[2*MK]; int ns=0; for(int k=0;k<MK;k++) if(active_(k)){ s[ns++]=P[k]; s[ns++]=Q[k]; } return paar(s,ns,randomize,NULL); }
static double total(int *gates){ int na=nact(); int g=slpcost(0)+2*na-1; if(gates) *gates=g; return g+LAM*unc; }
static void dump(const char *fn){
  FILE *f=fopen(fn,"w"); for(int k=0;k<MK;k++) if(active_(k)) fprintf(f,"%016llx %016llx\n",(unsigned long long)P[k],(unsigned long long)Q[k]); fclose(f);
}
int main(int argc,char**argv){
  if(argc<7){fprintf(stderr,"psa init(terms file or -) out seed iters T0 lam\n");return 1;}
  for(int v=0;v<32;v++){ COLS[v]=0; for(int u=0;u<32;u++) if(ROWS[u]>>v&1) COLS[v]|=1u<<u; }
  rs^=(u64)atoll(argv[3])*0x9E3779B97F4A7C15ULL; for(int i=0;i<20;i++) rnd();
  long iters=atol(argv[4]); double T0=atof(argv[5]); LAM=atof(argv[6]);
  unc=0; for(int u=0;u<32;u++) for(int v=0;v<32;v++){ cnt[u][v]=0; if(ROWS[u]>>v&1) unc++; }
  for(int k=0;k<MK;k++){P[k]=Q[k]=0;}
  if(strcmp(argv[1],"-")){ FILE *f=fopen(argv[1],"r"); int k=0; unsigned long long a,b; while(k<MK && fscanf(f,"%llx %llx",&a,&b)==2){ P[k]=a; Q[k]=b; if(!valid(P[k],Q[k])){fprintf(stderr,"invalid term %d\n",k);return 1;} term(k,1); k++; } fclose(f); }
  else { for(int v=0;v<32;v++){ P[v]=COLS[v]; Q[v]=1ULL<<(32+v); term(v,1);} }
  int g; double cost=total(&g); double best=1e9; int bestg=0;
  fprintf(stderr,"start gates %d unc %d\n",g,unc);
  for(long it=0;it<iters;it++){
    double T=T0*(1.0-(double)it/iters)+0.02;
    int k=rnd()%MK; u64 op=P[k], oq=Q[k]; u64 np=op, nq=oq;
    int mv=rnd()%10;
    if(mv<4){ // add literal
      int lit=rnd()%64; if(rnd()&1) np|=1ULL<<lit; else nq|=1ULL<<lit; if(np==op&&nq==oq) continue;
      if(!op||!oq){ // empty term: seed new simple rectangle cell
        int u=rnd()%32,v=rnd()%32; if(!(ROWS[u]>>v&1)) continue; np=1ULL<<u; nq=1ULL<<(32+v); }
    } else if(mv<8){ // remove literal
      if(!op&&!oq) continue; u64 s= (rnd()&1)? op:oq; if(!s) continue; int c=__builtin_popcountll(s); int t=rnd()%c; u64 bit=0; for(int b=0;b<64;b++) if(s>>b&1){ if(t--==0){bit=1ULL<<b;break;} }
      if(s==op) np&=~bit; else nq&=~bit;
    } else if(mv==8){ // move literal from term k to another term j
      int j=rnd()%MK; if(j==k||!active_(j)||!active_(k)) continue; int lit=rnd()%64; u64 bit=1ULL<<lit;
      int side=rnd()&1; u64 *src= side? &P[k]:&Q[k]; if(!(*src&bit)) continue;
      u64 jp=P[j], jq=Q[j]; int side2=rnd()&1; u64 njp=jp,njq=jq; if(side2) njp|=bit; else njq|=bit; if(njp==jp&&njq==jq) continue;
      if(!valid(njp,njq)) continue;
      u64 nk_p= side? (op&~bit):op, nk_q= side? oq:(oq&~bit);
      term(k,-1); term(j,-1); P[k]=nk_p; Q[k]=nk_q; P[j]=njp; Q[j]=njq; term(k,1); term(j,1);
      int g2; double nc=total(&g2); double d=nc-cost;
      if(d<=0 || rndf()<exp(-d/T)){ cost=nc; g=g2; }
      else { term(k,-1); term(j,-1); P[k]=op; Q[k]=oq; P[j]=jp; Q[j]=jq; term(k,1); term(j,1); }
      goto rec;
    } else { // swap sides of one literal within term (P->Q)
      int lit=rnd()%64; u64 bit=1ULL<<lit; if(op&bit){ np&=~bit; nq|=bit; } else if(oq&bit){ nq&=~bit; np|=bit; } else continue;
    }
    if(np&&nq&&!valid(np,nq)) continue;
    term(k,-1); P[k]=np; Q[k]=nq; term(k,1);
    { int g2; double nc=total(&g2); double d=nc-cost;
      if(d<=0 || rndf()<exp(-d/T)){ cost=nc; g=g2; }
      else { term(k,-1); P[k]=op; Q[k]=oq; term(k,1); } }
  rec:
    if(unc==0 && g<best){ best=g; bestg=g; dump(argv[2]); fprintf(stderr,"[%ld] best %d (K=%d)\n",it,g,nact()); }
  }
  printf("%d\n",bestg); return 0;
}
'''
CGP_C = r'''// CGP (1+lambda) gate-count minimisation for the promise-represent problem.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#define W 16
#define NI 64
#define MAXG 1024
typedef uint64_t u64;
static u64 TGT[W];
static u64 IN[NI][W];
typedef struct { int N; unsigned char f[MAXG]; short a[MAXG], b[MAXG]; int out; } G; // node id = NI + i
static u64 rs=88172645463325252ULL;
static inline u64 rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static u64 VAL[NI+MAXG][W];
static unsigned char act[NI+MAXG];
static int active(G *g){ // marks active, returns count of active gates
  memset(act,0,NI+g->N);
  act[g->out]=1; int c=0;
  for(int i=g->N-1;i>=0;i--){ int id=NI+i; if(!act[id]) continue; c++; act[g->a[i]]=1; if(g->f[i]!=3) act[g->b[i]]=1; }
  return c;
}
static int evalg(G *g){ // returns errors; assumes act computed
  for(int i=0;i<g->N;i++){ int id=NI+i; if(!act[id]) continue; u64 *A=VAL[g->a[i]], *B=VAL[g->b[i]], *R=VAL[id];
    switch(g->f[i]){ case 1: for(int w=0;w<W;w++) R[w]=A[w]&B[w]; break; case 2: for(int w=0;w<W;w++) R[w]=A[w]|B[w]; break; default: for(int w=0;w<W;w++) R[w]=~A[w]; }
  }
  int e=0; for(int w=0;w<W;w++) e+=__builtin_popcountll(VAL[g->out][w]^TGT[w]); return e;
}
static void save(const char *fn,G *g){
  active(g); static int map[NI+MAXG]; for(int i=0;i<NI;i++) map[i]=i; int n=0;
  for(int i=0;i<g->N;i++) if(act[NI+i]) map[NI+i]=NI+n++;
  char tmp[600]; snprintf(tmp,600,"%s.tmp",fn); FILE *f=fopen(tmp,"w"); fprintf(f,"%d %d\n",n,map[g->out]);
  for(int i=0;i<g->N;i++) if(act[NI+i]){ if(g->f[i]==3) fprintf(f,"not %d -1\n",map[g->a[i]]); else fprintf(f,"%s %d %d\n",g->f[i]==1?"and":"or",map[g->a[i]],map[g->b[i]]); }
  fclose(f); rename(tmp,fn);
}
static void randnode(G *g,int i){ int lim=NI+i; g->f[i]=1+rnd()%3; if(rnd()%4) g->f[i]= 1+rnd()%2; g->a[i]=rnd()%lim; g->b[i]=rnd()%lim; }
int main(int argc,char**argv){
  if(argc<6){fprintf(stderr,"cgp init out seed seconds slackfactor [lambda]\n");return 1;}
  unsigned ROWS[32]={0xf2bc7573, 0x59092a9d, 0xc443c8c2, 0x1d58483a, 0x1b7982bb, 0xa8882290, 0x74d0fd25, 0x00041402, 0x6879a8a5, 0xd65a6f1a, 0x7455fba3, 0x24508c1a, 0x3050a010, 0x68658ad5, 0x10a00735, 0x5e55690e, 0x6ab804a4, 0x00b01458, 0xf447af40, 0x00b07351, 0xa88c2290, 0x40d9d7bc, 0x40241486, 0x54058082, 0x6891afb9, 0x50ddf79e, 0x553080be, 0x28c180f8, 0x6ac8ae9e, 0xa6860e80, 0x20508c18, 0x10f01050};
  for(int w=0;w<W;w++) TGT[w]=0;
  for(int u=0;u<32;u++) TGT[u>>1] |= ((u64)ROWS[u])<<((u&1)*32);
  for(int u=0;u<32;u++){ for(int w=0;w<W;w++) IN[u][w]=0; IN[u][u>>1]=(u&1)?0xFFFFFFFF00000000ULL:0xFFFFFFFFULL; }
  for(int v=0;v<32;v++) for(int w=0;w<W;w++) IN[32+v][w]=(1ULL<<v)|(1ULL<<(v+32));
  for(int i=0;i<NI;i++) memcpy(VAL[i],IN[i],sizeof(IN[i]));
  rs^=(u64)atoll(argv[3])*0x9E3779B97F4A7C15ULL; for(int i=0;i<20;i++) rnd();
  double secs=atof(argv[4]); double slack=atof(argv[5]); int lambda= argc>6? atoi(argv[6]):4;
  // load
  FILE *f=fopen(argv[1],"r"); int n,o; fscanf(f,"%d %d",&n,&o);
  int fop[MAXG],fa[MAXG],fb[MAXG];
  for(int i=0;i<n;i++){ char s[16]; int a,b=-1; fscanf(f,"%15s %d %d",s,&a,&b); fop[i]= s[0]=='a'?1: s[0]=='o'?2:3; fa[i]=a; fb[i]=b; }
  fclose(f);
  int N=(int)(n*slack); if(N>MAXG) N=MAXG; if(N<n) N=n;
  // place original gates at evenly spread positions
  G p; p.N=N; static int posmap[NI+MAXG]; for(int i=0;i<NI;i++) posmap[i]=i;
  int slot[MAXG]; for(int i=0;i<n;i++) slot[i]= (int)((long)i*N/n);
  for(int i=0;i<N;i++) p.f[i]=0;
  for(int i=0;i<n;i++) posmap[NI+i]=NI+slot[i];
  for(int i=0;i<n;i++){ int s=slot[i]; p.f[s]=fop[i]; p.a[s]=posmap[fa[i]]; p.b[s]= fop[i]==3? posmap[fa[i]] : posmap[fb[i]]; }
  for(int i=0;i<N;i++) if(!p.f[i]) randnode(&p,i);
  p.out=posmap[o];
  int ps=active(&p); int e=evalg(&p);
  if(e){ fprintf(stderr,"init wrong %d\n",e); return 1; }
  fprintf(stderr,"start %d (N=%d)\n",ps,N);
  int best=ps; save(argv[2],&p);
  time_t t0=time(NULL); long gen=0;
  G c, bc;
  while(1){
    gen++;
    if((gen&1023)==0 && difftime(time(NULL),t0)>secs) break;
    int bsz=1<<30; int found=0;
    for(int l=0;l<lambda;l++){
      c=p;
      // single active mutation
      active(&c);
      for(int tries=0;tries<1000;tries++){
        int r=rnd()%(3*N+1);
        if(r==3*N){ int no=rnd()%(NI+N); if(no!=c.out){ c.out=no; break;} continue; }
        int i=r/3, gene=r%3; int lim=NI+i;
        if(gene==0){ int nf=1+rnd()%3; if(nf==c.f[i]) continue; c.f[i]=nf; }
        else if(gene==1){ int na=rnd()%lim; if(na==c.a[i]) continue; c.a[i]=na; }
        else { int nb=rnd()%lim; if(nb==c.b[i]) continue; c.b[i]=nb; if(c.f[i]==3) { continue; } }
        if(act[NI+i]) break;
      }
      int sz=active(&c);
      if(sz>ps) continue;
      int er=evalg(&c);
      if(er) continue;
      if(sz<bsz || (sz==bsz && (rnd()&1))){ bsz=sz; bc=c; found=1; }
    }
    if(found){ p=bc; ps=bsz; if(ps<best){ best=ps; save(argv[2],&p); fprintf(stderr,"[%ld] best %d\n",gen,best); fflush(stderr);} }
  }
  save(argv[2],&p); printf("%d\n",best); return 0;
}
'''
SEED_TERMS = '''0080040400100221 8008201012204040
0010041800004441 1000012141020800
1100000802200441 401051000400001a
2040a00002400081 0004140210001040
0041804020010001 2280040008042004
48000804800a0001 00a010401000a002
1405000012200022 0008029000848000
0410001201000502 0001200104000808
000080a80000044c 044048000d200000
2880000020040204 8402080009100020
0818002000800404 4000808200010118
2004028082001200 1050001020102400
4020008410008210 0240000a04410100
2040880000086010 0020021140040900
2000088002088200 0010610001000000
0041809801040200 400027000a200010
'''
SEED_CIRC = '''190 253
or 9 39
or 8 33
or 6 10
or 47 54
or 35 66
or 1 52
or 4 62
or 24 58
or 25 60
or 5 20
or 11 32
or 28 53
or 0 34
or 0 67
or 2 18
or 2 37
or 3 26
or 3 47
or 4 21
or 10 52
or 12 50
or 13 36
or 13 61
or 14 42
or 15 75
or 16 65
or 19 43
or 21 72
or 23 51
or 25 45
or 28 73
or 30 74
or 36 41
or 36 68
or 39 48
or 40 42
or 40 46
or 43 71
or 44 69
or 61 64
or 18 96
or 27 101
or 55 78
or 9 55
or 107 73
or 108 76
or 109 99
or 6 14
or 111 21
or 112 51
or 113 63
or 114 93
or 15 19
or 45 100
or 0 52
or 118 87
or 119 97
or 17 24
or 121 37
or 122 72
or 123 95
or 0 56
or 125 68
or 126 91
or 70 80
or 128 100
or 129 102
or 7 22
or 131 61
or 132 77
or 133 93
or 28 42
or 135 44
or 136 65
or 137 66
or 138 84
or 16 29
or 140 38
or 141 48
or 142 77
or 106 27
or 144 42
or 145 57
or 146 86
or 17 31
or 148 59
or 149 62
or 150 76
or 151 90
or 38 55
or 153 85
or 154 88
or 155 102
or 1 48
or 157 50
or 158 58
or 159 91
or 160 94
or 104 15
or 162 92
or 65 69
or 164 71
or 15 43
or 166 45
or 167 74
or 168 80
or 169 98
or 64 79
or 171 81
or 172 97
or 105 26
or 174 46
or 175 82
or 48 54
or 106 29
or 178 59
or 179 103
or 105 49
or 181 63
or 182 94
or 25 59
or 184 79
or 185 83
or 186 92
or 70 81
or 188 89
or 189 98
or 3 25
or 191 31
or 192 41
or 193 84
or 194 103
or 20 29
or 196 54
or 197 60
or 198 83
or 199 85
or 34 64
or 201 70
or 202 88
or 22 26
or 204 35
or 205 54
or 206 57
or 207 89
or 4 67
or 209 86
or 210 87
or 211 90
or 104 75
or 213 95
or 9 18
or 215 24
or 216 25
or 41 45
or 218 62
or 219 99
and 110 115
and 120 124
and 127 130
and 134 139
and 143 147
and 152 156
and 161 163
and 165 170
and 173 176
and 180 183
and 187 190
and 195 200
and 203 208
and 212 214
and 116 117
and 82 177
and 217 220
or 221 222
or 223 224
or 225 226
or 227 228
or 229 230
or 231 232
or 233 234
or 235 236
or 238 239
or 240 241
or 242 243
or 244 245
or 246 247
or 248 249
or 250 251
or 252 237
'''
TARGET = 140

def verify(gates, out):
    """gates: list of (op, [inputs]); exact Lean semantics on all 1024 one-a-side inputs. Returns #errors."""
    for i, (op, ins) in enumerate(gates):
        if op not in ('and', 'or', 'not') or len(ins) > 2 or any(j >= 64 + i or j < 0 for j in ins):
            return 10**9
    bad = 0
    for u in range(32):
        for v in range(32):
            w = [i == u for i in range(32)] + [j == v for j in range(32)]
            for op, ins in gates:
                vals = [w[i] if i < len(w) else False for i in ins]
                w.append(all(vals) if op == 'and' else any(vals) if op == 'or' else not (vals[0] if vals else False))
            if (w[out] if out < len(w) else False) != bool((ROWS[u] >> v) & 1):
                bad += 1
    return bad

def paar_circuit(terms, rng):
    """terms: list of (P,Q) 64-bit masks. Build OR chains with randomized Paar, then ANDs and OR tree."""
    sets = []
    for p, q in terms:
        sets += [p, q]
    uniq = list(dict.fromkeys(sets))
    reps = {s: set(b for b in range(64) if s >> b & 1) for s in uniq}
    gates = []
    while True:
        cnt = {}
        for s, r in reps.items():
            if len(r) < 2: continue
            l = sorted(r)
            for i in range(len(l)):
                for j in range(i + 1, len(l)):
                    cnt[(l[i], l[j])] = cnt.get((l[i], l[j]), 0) + 1
        if not cnt: break
        m = max(cnt.values()); a, b = rng.choice([k for k, c in cnt.items() if c == m])
        nid = 64 + len(gates); gates.append(('or', [a, b]))
        for r in reps.values():
            if a in r and b in r: r.discard(a); r.discard(b); r.add(nid)
    wire = {s: next(iter(r)) for s, r in reps.items()}
    tl = []
    for p, q in terms:
        nid = 64 + len(gates); gates.append(('and', [wire[p], wire[q]])); tl.append(nid)
    while len(tl) > 1:
        a = tl.pop(0); b = tl.pop(0); nid = 64 + len(gates); gates.append(('or', [a, b])); tl.append(nid)
    return gates, tl[0]

def load_circ_text(txt):
    L = txt.strip().split('\n'); n, o = map(int, L[0].split()); g = []
    for l in L[1:1 + n]:
        p = l.split(); g.append((p[0], [int(x) for x in p[1:] if int(x) >= 0]))
    return g, o

def circ_text(g, o):
    return '%d %d\n' % (len(g), o) + ''.join('%s %d %d\n' % (op, ins[0], ins[1] if len(ins) > 1 else -1) for op, ins in g)

class Best:
    def __init__(self, path):
        self.path = path; self.size = 10**9; self.g = None; self.o = None
        if os.path.exists(path):
            try:
                d = json.load(open(path)); g = [(x[0], x[1]) for x in d['gates']]
                if verify(g, d['output']) == 0: self.size, self.g, self.o = len(g), g, d['output']
            except Exception: pass
    def offer(self, g, o, how):
        if len(g) >= self.size: return False
        if verify(g, o) != 0: return False
        self.size, self.g, self.o = len(g), g, o
        tmp = self.path + '.tmp'
        json.dump({'step': 'L-meta-20581', 'gate_count': len(g), 'output': o, 'gates': [[op, ins] for op, ins in g],
                   'found_by': how, 'time': time.strftime('%Y-%m-%d %H:%M:%S'),
                   'semantics': 'wires 0..31 = x_u (left), 32..63 = y_v (right), gate j -> wire 64+j'}, open(tmp, 'w'))
        os.replace(tmp, self.path)
        print('BEST %d (%s)' % (len(g), how), flush=True)
        return True

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--hours', type=float, default=5.5)
    ap.add_argument('--jobs', type=int, default=max(1, (os.cpu_count() or 2)))
    ap.add_argument('--out', default=os.path.join(HERE, 'results', '20581_best.json'))
    a = ap.parse_args(); deadline = time.time() + a.hours * 3600
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    wd = tempfile.mkdtemp(prefix='s20581_')
    cc = shutil.which('gcc') or shutil.which('cc') or shutil.which('clang')
    for name, src in (('psa', PSA_C), ('cgp', CGP_C)):
        open(os.path.join(wd, name + '.c'), 'w').write(src)
        subprocess.run([cc, '-O3', '-march=native', '-o', os.path.join(wd, name), os.path.join(wd, name + '.c'), '-lm'], check=True)
    best = Best(a.out); rng = random.Random()
    g, o = load_circ_text(SEED_CIRC); best.offer(g, o, 'fleet L_meta_19143 circuit')
    pool = [[int(x, 16) for x in l.split()] for l in SEED_TERMS.strip().split('\n')]
    pool_terms = [(pool, 0)]
    def terms_to_best(terms, how):
        for _ in range(30):
            g, o = paar_circuit(terms, rng); best.offer(g, o, how)
    terms_to_best([tuple(t) for t in pool], 'seed terms')
    tf = os.path.join(wd, 'bestterms.txt')
    open(tf, 'w').write('\n'.join('%016x %016x' % tuple(t) for t in pool) + '\n')
    rnd = 0
    while time.time() < deadline and best.size > TARGET:
        rnd += 1; procs = []
        njobs = a.jobs
        for j in range(njobs):
            seed = rng.randrange(1, 10**9)
            if j < max(1, njobs // 3):  # CGP polish of current best circuit
                cf = os.path.join(wd, 'cin%d.txt' % j); open(cf, 'w').write(circ_text(best.g, best.o))
                outf = os.path.join(wd, 'cout%d.txt' % j)
                secs = int(max(30, min(1200, deadline - time.time() - 60)))
                procs.append(('cgp', outf, subprocess.Popen([os.path.join(wd, 'cgp'), cf, outf, str(seed), str(secs), '1.5', '4'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)))
            else:
                outf = os.path.join(wd, 'pout%d.txt' % j)
                init = tf if rng.random() < 0.8 else '-'
                T0 = rng.choice([0.3, 0.5, 0.6, 0.8, 1.0]); lam = rng.choice([2.5, 3.0, 3.5])
                its = rng.choice([3000000, 6000000, 12000000])
                procs.append(('psa', outf, subprocess.Popen([os.path.join(wd, 'psa'), init, outf, str(seed), str(its), str(T0), str(lam)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)))
        for kind, outf, p in procs:
            try:
                p.wait(timeout=max(1, deadline - time.time()))
            except subprocess.TimeoutExpired:
                p.kill()
            if not os.path.exists(outf): continue
            try:
                if kind == 'cgp':
                    g, o = load_circ_text(open(outf).read()); best.offer(g, o, 'cgp polish')
                else:
                    terms = [tuple(int(x, 16) for x in l.split()) for l in open(outf).read().strip().split('\n') if l.strip()]
                    before = best.size
                    terms_to_best(terms, 'cross-pair SA + Paar')
                    if best.size < before:
                        open(tf, 'w').write('\n'.join('%016x %016x' % t for t in terms) + '\n')
            except Exception as e:
                print('warn', e, flush=True)
        print('round %d best %d' % (rnd, best.size), flush=True)
    print('DONE best %d' % best.size, flush=True)

if __name__ == '__main__':
    main()
