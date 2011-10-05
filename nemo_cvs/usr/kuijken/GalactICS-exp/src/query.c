/*
 * This file contains four functions for prompting single value input, 
 * one for each of the data types char, int, float and double.
 * The functions are cquery, iquery, fquery and dquery.
 * Each function has two arguments: a prompt and a pointer variable of the
 * appropriate data type.  
 * The function will display the prompt (e.g. "Enter x") and a default
 * value.  The user can either enter a new value or accept the default
 * by hitting a return.
 * 
 * Example
 * main()
 * {
 *		char a[50];
 *		float x=5;
 *		
 *		strcpy(a,"hello there");
 *		cquery("Enter a phrase",a);
 *		fquery("Enter the value to be squared",&x);
 *      fprintf(stdout,"%s %g squared = %g\n",a,x*x);
 * }
 */

#include <stdio.h>
double atof();

cquery(prompt,c)
char *c;
char *prompt;
{
	char *p, line[100], newprompt[200];
	strcpy(newprompt,prompt);
	strcat(newprompt," [%s]: ");
	fprintf(stderr,newprompt,c);
	fgets(line,80,stdin);
	if(*line != '\n')
	{
		for(p=line; *p != '\n' && *p != '#'; p++);
		if( *p == '#' ) {
			p--;
			while( *p == ' ' || *p == '\t' ) p--;
			p++;
		}
		*p = '\0';
		
		strcpy(c,line);
	}
}

dquery(prompt,d)
double *d;
char prompt[100];
{
	char *p, line[100], newprompt[200];
	strcpy(newprompt,prompt);
	strcat(newprompt," [%g]: ");
	fprintf(stderr,newprompt,*d);
	fgets(line,80,stdin);
	if(*line != '\n')
	{
		for(p=line; *p !='\n' && *p !='#'; p++);
			*p = '\0';
		*d = atof(line);
	}
}

fquery(prompt,f)
float *f;
char prompt[100];
{
	char *p, line[100], newprompt[200];
	strcpy(newprompt,prompt);
	strcat(newprompt," [%g]: ");
	fprintf(stderr,newprompt,*f);
	fgets(line,80,stdin);
	if(*line != '\n')
	{
		for(p=line; *p!='\n' && *p!='#'; p++);
			*p = '\0';
		*f = (float) atof(line);
	}
}

iquery(prompt,i)
int *i;
char prompt[100];
{
	char *p, line[100], newprompt[200];
	strcpy(newprompt,prompt);
	strcat(newprompt," [%d]: ");
	fprintf(stderr,newprompt,*i);
	fgets(line,80,stdin);
	if(*line != '\n')
	{
		for(p=line; *p !='\n' && *p !='#'; p++);
			*p = '\0';
		*i = atoi(line);
	}
}
