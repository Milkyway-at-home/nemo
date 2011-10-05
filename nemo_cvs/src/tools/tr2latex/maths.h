/*
** tr2latex - troff to LaTeX converter
** $Id: maths.h,v 1.1.1.1 1992/04/27 15:20:59 teuben Exp $
** COPYRIGHT (C) 1987 Kamal Al-Yahya, 1991,1992 Christian Engel
**
** Module: maths.h
**
** This file contains a list of the words that have simple
** correspondence in the two languages.
** Do not put here words that require action (like sub).
** If the word is identical in the two languages (except for TeX's backslash),
** put it in simil.h
*/

struct math_equiv {
	char *troff_symb, *tex_symb;
} math[] = {
/*	troff name		TeX name		*/

	"~",			"\\ ",
	"^",			"\\,",
	"above",		"\\cr",
	"ccol",			"\\matrix",
	"cpile",		"\\matrix",
	"fat",			"",
	"grad",			"\\nabla",
	"half",			"{1\\over 2}",
	"inf",			"\\infty",
	"inter",		"\\cap",
	"lcol",			"\\matrix",
	"lineup",		"",
	"lpile",		"\\matrix",
	"mark",			"",
	"nothing",		"",
	"pile",			"\\matrix",
	"rcol",			"\\matrix",
	"rpile",		"\\matrix",
	"union",		"\\cup"
};
