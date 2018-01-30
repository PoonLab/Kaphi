%{
#include <stdio.h>

#include "igraph.h"

#include "tree.h"
#include "newick_parser.h"

#pragma GCC diagnostic ignored "-Wswitch-default"
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wredundant-decls"
%}

%%
[\n\t\r ] ;
#[^\n]*\n   ;
\(  return LPAREN;
\)  return RPAREN;
:   return COLON;
;   return SEMICOLON;
,   return COMMA;
[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?  yylval.number = atof(yytext); return NUMBER;
[^:;(),]+   yylval.string = yytext; return STRING;
%%
