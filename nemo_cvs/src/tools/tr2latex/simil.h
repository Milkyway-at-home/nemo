/*
** tr2latex - troff to LaTeX converter
** $Id: simil.h,v 1.1.1.1 1992/04/27 15:20:56 teuben Exp $
** COPYRIGHT (C) 1987 Kamal Al-Yahya, 1991,1992 Christian Engel
**
** Module: simil.h
**
** This file contains a list of math words that are similar in the
** two languages (in fact identical except for TeX's backslah).
** If I overlooked anything out, it can be put here
** Do NOT put here words that are similar but require some action (like over)
*/

char *sim_list[] =
{
"alpha",    "approx",    "beta",     "cdot",     "chi",      "cos",
"cosh",     "cot",       "coth",     "delta",    "epsilon",  "eta",
"exp",      "gamma",     "int",      "kappa",    "lambda",   "lim",
"log",      "matrix",    "mu",       "nu",       "omega",    "partial",
"phi",      "pi",        "prime",    "prod",     "psi",      "rho",
"sigma",    "sin",       "sinh",     "sqrt",     "sum",      "tan",
"tanh",     "tau",       "theta",    "times",    "xi",       "zeta"
};
