/***************************************************************/
/* File: loadobj.c                                             */
/* Last modified on Wed Aug 20 09:22:36 1986 by roberts        */
/* Last modified on Sat Dec 06          1986 by josh           */
/* Last modified on Wed May 17 03:17:00 1989 by Peter (COFF)   */
/*		        July 1990 - some hacking  by PJT       */
/*                  November 1991 - AIX stuff                  */
/*	AIX version with XCOFF cloned off 3b1 version          */
/* % make loadobjtest MACH=AIX EL=-lld                         */
/*	and kept in pace with 3b1 using #if defined(aix)       */
/* Externals called from the libld.a library:                  */
/*   ldopen  ldnshread  ldnsseek  ldnrseek  ldohseek  ldtbseek */
/* On AIX there is a difference (due to alignment) of sizeof   */
/* and the actual size of the object, eg RELSZ vs. sizeof item */
/***************************************************************/

#include <stdinc.h>
#include <a.out.h>	/* also includes things as filehdr.h etc... in SYS5 */
#include <ldfcn.h>      /* SYS5 */
#include <strlib.h>
#include <filefn.h>
#include <loadobj.h>

/* next one is for silly 3b1: see a.out.h -> syms.h and nlist.h */
/* #undef n_name */

/***************************************************************/
/* Types for binary tree symbol table                          */
/***************************************************************/

typedef struct _ste {
   struct nlist stdata;
   struct _ste *before,*after;
} *symtree;

typedef struct nlist *symbol;

/***************************************************************/
/* Package variables                                           */
/*                                                             */
/*     The use of these variables describes their use to a     */
/* reasonably sufficient level of detail.  The only subtlety   */
/* is that localtable is used as a dynamically-allocated       */
/* array of symbol blocks.                                     */
/***************************************************************/

static symtree symbase;

#if 0
static stream infile;				/* not used in SYS5 */
static struct exec header;			/* not used in SYS5 */
#endif
static LDFILE *ldptr;             /* replaces infile/header in SYS5 */
static FILHDR   header;                          /* this is in SYS5 */
static AOUTHDR  oheader;                         /* this is in SYS5 */

static char *stringtable;    /* first 4 bytes = length ; then strings */
static struct nlist *localtable;
static long nsymbols;

static char *textsegment;	/* points to text of object file */
static char *datasegment;	/*	     data                */
static char *bsssegment;	/*	     bss		 */

/***************************************************************/
/* Local function declarations                                 */
/***************************************************************/

static symbol lookup(/* name */);
static symbol enter(/* sym */);
static symbol findentry(/* name, eptr, value */);
static void readstrings();
static void emptystrings();
static void readsymbols(/* reloc */);
static void processrelocation(/* size, segment */);
static string savename( /* char[8] */);



/***************************************************************/
/* loadobj(pathname);                                          */
/*                                                             */
/*     Loads the object file indicated by pathname.  If        */
/* successful, nothing is returned.  Errors are handled        */
/* by internal calls on the error routine.  The global         */
/* symbols from this file are added to the runtime database.   */
/***************************************************************/

void loadobj(pathname)
string pathname;
{
    char *newname;
    SCNHDR texthead, datahead, bsshead;

    newname = pathfind("", pathname);
    if (newname == NULL)
        error("loadobj: cannot find %s\n",pathname);
    ldptr = ldopen(newname,NULL);
    if (ldptr == NULL)
	error("loadobj: cannot ldopen object file %s", newname);
    dprintf(4,"\nloadobj : object file = %s\n",newname);
    dprintf(5,"              type = 0x%x\n",TYPE(ldptr));
    dprintf(5,"            offset = 0x%x\n",OFFSET(ldptr));
    dprintf(5,"         #sections = %d\n",HEADER(ldptr).f_nscns);
    dprintf(5,"         #symbols  = %d\n",HEADER(ldptr).f_nsyms);
    dprintf(5,"  aout opthdr size = %d\n",HEADER(ldptr).f_opthdr);

    if (OFFSET(ldptr) != 0)     /* don't allow archives files for now */
        error("loadobj: COFF nonzero offset for filehdr; .a file?");
    if (HEADER(ldptr).f_opthdr != 0)
        warning("loadobj: this COFF has an optional Unix System header");

    ldnshread (ldptr, ".text", &texthead);
    ldnshread (ldptr, ".data", &datahead);
    ldnshread (ldptr, ".bss",   &bsshead);
    dprintf(5,"   text -> 0x%x dat@ 0x%x(0x%x) rel@ 0x%x(%d) \n",
		texthead.s_paddr, texthead.s_scnptr,texthead.s_size,
		texthead.s_relptr,texthead.s_nreloc);
    dprintf(5,"   data -> 0x%x dat@ 0x%x(0x%x) rel@ 0x%x(%d) \n",
		datahead.s_paddr, datahead.s_scnptr,datahead.s_size,
		datahead.s_relptr,datahead.s_nreloc);
    dprintf(5,"    bss -> 0x%x dat@ 0x%x(0x%x) rel@ 0x%x(%d)\n", 
		bsshead.s_paddr,  bsshead.s_scnptr,bsshead.s_size,
		bsshead.s_relptr, bsshead.s_nreloc);
    textsegment = 
       getmem( (int) texthead.s_size + datahead.s_size + bsshead.s_size);
    datasegment = textsegment + texthead.s_size;
    bsssegment = datasegment + datahead.s_size;
    dprintf(5,"   Textsegment allocated at 0x%x\n",textsegment);
    if (texthead.s_size > 0) {           /* there may be no text... */
    	if (ldnsseek(ldptr, ".text") == FAILURE) 
        	error("loadobj: cannot ldnsseek to .text\n");
    	FREAD(textsegment, 1, texthead.s_size,ldptr);
    }
    if (datahead.s_size > 0) {           /* there may be no data ... */
        if (ldnsseek(ldptr, ".data") == FAILURE)
            error("loadobj: cannot ldnsseek to .data\n");
        FREAD(datasegment, 1, datahead.s_size,ldptr);
    }
    readstrings();
    readsymbols(TRUE);

    ldnrseek(ldptr,".text");                        /* seek to text part */
    processrelocation(texthead.s_nreloc, textsegment);/* and process it */
    ldnrseek(ldptr,".data");                        /* seek to data part */
    processrelocation(datahead.s_nreloc, datasegment);/* and process it */

    free((char *) localtable);                      /* free locally used -- */
    free((char *) stringtable);                     /* -- tables */
}



/***************************************************************/
/* fn = findfn(fnname);                                        */
/* (*fn)();                                                    */
/*                                                             */
/*      The findfn routine looks up fnname in the symbol table */
/* and returns a function pointer by which that function may   */
/* be called.  If the name is not defined, NULL is returned.   */
/*      PJT comments:                                          */ 
/*      Shouldn't findfn() know about the naminbg convention   */
/*      of symbols of your host system (underscores, case)     */
/***************************************************************/

proc findfn(fnname)
string fnname;
{
    symbol sym;

    sym = lookup(fnname);
    return ((sym == NULL) ? NULL : (proc) sym->n_value);
}



/***************************************************************/
/* mysymbols(progname);                                        */
/*                                                             */
/*     Loads the symbols defined in progname, which is usually */
/* argv[0] from the main program.  In case argv[0] is not the  */
/* complete name, the path environment variable is used.       */
/***************************************************************/

void mysymbols(progname)
string progname;
{
    char *newname;

    newname = pathfind(getenv("PATH"), progname);
    if (newname == NULL)
        error("mysymbols: can't find %s along PATH",progname);
    ldptr = ldopen(newname,NULL);
    if (ldptr == NULL)
	error("mysymbols: cant ldopen executable %s", newname);
    dprintf(4,"\nmysymbols : exec file = %s\n",newname);
    dprintf(5,"              type = 0x%x\n",TYPE(ldptr));
    dprintf(5,"            offset = 0x%x\n",OFFSET(ldptr));
    dprintf(5,"         #sections = %d\n",HEADER(ldptr).f_nscns);
    dprintf(5,"         #symbols  = %d\n",HEADER(ldptr).f_nsyms);
    dprintf(5,"  aout opthdr size = %d\n",HEADER(ldptr).f_opthdr);

    if (OFFSET(ldptr) != 0)
        error("mysymbols: COFF nonzero offset for filehdr??? \n");
    if (HEADER(ldptr).f_opthdr == 0)
        error("mysymbols: COFF has no Unix System header??\n");
    if (ldohseek(ldptr) == FAILURE)
        error("mysymbols: cannot ldohseek to Unix System header\n");
    FREAD((char *) &oheader, sizeof(AOUTHDR), 1, ldptr);
    dprintf(5,"     AOUTHDR :  magic = 0x%x\n",oheader.magic);
    dprintf(5,"             textsize = %d\n",oheader.tsize);
    dprintf(5,"             datasize = %d\n",oheader.dsize);
    dprintf(5,"              bsssize = %d\n",oheader.bsize);
    dprintf(5,"            entry pt = 0x%x\n",oheader.entry);
    dprintf(5,"            text @ 0x%x\n",oheader.text_start);
    dprintf(5,"            data @ 0x%x\n",oheader.data_start);
    dprintf(5,"            o_toc @ 0x%x\n",oheader.o_toc);

    readstrings();
    readsymbols(FALSE);
    free(stringtable);
}



/***************************************************************/
/* readstrings();                                              */
/*                                                             */
/*     Reads in the complete string table from the a.out file. */
/* This storage is freed at the end of the loadfn call.        */
/*                          (mysymbols or loadobj)             */
/***************************************************************/

static void readstrings()
{
    int size, err;

    dprintf(5,"readstrings:\n");
    if (SYMESZ != AUXESZ)
        error("readstrings: how can I find strings????? check COFF code\n");
    if (ldtbseek(ldptr) == FAILURE) { /* I'm not sure if we should allow this */
					/* maybe just error() out here */
        dprintf(5,"   warning: (COFF) cannot even seek to symbol section\n");
        emptystrings();
        return;
    }
    dprintf(5,"   seeking to symbols    ftell: @ 0x%x\n",FTELL(ldptr));
    err = FSEEK(ldptr, SYMESZ * HEADER(ldptr).f_nsyms, CURRENT);
    dprintf(5,"   skipping symbols, now ftell: @ 0x%x fseek-> %d\n",
				FTELL(ldptr),err);
    if (err != 0) {			/* expected fseek to work */
        emptystrings();
        return;
    }
    err=FREAD((char *) &size, sizeof size, 1, ldptr);
    if (err != 1) {			/* expected only 1 item to read */
        dprintf(5,"   no stringtable, fread returned %d\n",err);
        emptystrings();         /* assume empty */
        return;
    }
    dprintf(5,"   size of stringtable = %d  fread-> %d\n",size,err);
    stringtable = (string) getmem(size);
    dprintf(5,"   stringtable allocated @ 0x%x\n",stringtable);
    err=FSEEK(ldptr, -(sizeof size), CURRENT);      /* go back */
    dprintf(5,"   reading strings @ 0x%x fseek-> %d\n",FTELL(ldptr),err);
    err=FREAD(stringtable, 1, size, ldptr);
    dprintf(5,"   read stringtable, fread-> %d\n",err);
    if (err != size)
        error("readstrings: read incomplete stringtable, size,err=%d %d\n",
                size,err);
}

static void emptystrings()	/* declare empty stringtable */
{
    stringtable = (string) getmem(5);
    dprintf(5,"   emptystrings : assume dummy empty stringable anyhow\n");
    dprintf(5,"   stringtable declared @ 0x%x\n",stringtable);
    *(stringtable+0) = 0;     /* PJT: fails if sizeof(int) != 4 */
    *(stringtable+1) = 0;
    *(stringtable+2) = 0;
    *(stringtable+3) = 4;
    *(stringtable+4) = 0;     /* just in case zero terminate it */
}




/***************************************************************/
/* readsymbols(reloc);                                         */
/*                                                             */
/*     Reads in all of the symbols and defines the external    */
/* symbols in the symtab tree.  If reloc is TRUE, the value    */
/* of each symbol is relocated relative to the start of the    */
/* text segment, and all the symbols are stored in localtable  */
/* for relocation.  The localtable storage should is released  */
/* at the end of each loadfn.                                  */
/***************************************************************/

static void readsymbols(reloc)
bool reloc;
{
    struct syment symentry;         /* SYS5 temporary symbol entry */
    struct nlist entryblk;
    symbol entry;
    int i;
    static int numaux;

    nsymbols = HEADER(ldptr).f_nsyms;
    dprintf(5,"\nreadsymbols: nsymbols=%d\n",nsymbols);
    dprintf(5,"SYMESZ=%d  sizeof symentry=%d\n",
            SYMESZ, sizeof(symentry));
    if (reloc) localtable = (symbol) getmem((int) nsymbols*SYMESZ);
    ldtbseek(ldptr);                        /* get to symbol table */
    numaux = 0;
    for (i = 0; i < nsymbols; i++)  {
        entry = (reloc) ? &localtable[i] : &entryblk;
        dprintf(10,"  reading new symbol @ 0x%x\n",FTELL(ldptr));
/* 	FREAD((char *) &symentry, sizeof symentry, 1, ldptr);/* syment/auxent */
	FREAD((char *) &symentry, SYMESZ, 1, ldptr);/* syment/auxent */
        if (numaux > 0) {
/*            if(numaux != 1) warning("numaux = %d != 1 ",numaux);  */
/*            numaux = 0;       /* at most one numaux, so reset, no countdown */
            numaux--;
            if (reloc) entry->n_name = NULL;        /* flag bad, in case */
            continue;                         /* no further processing needed */
        } else
            numaux = symentry.n_numaux;
        if (symentry.n_zeroes == 0) 
            entry->n_name = stringtable + symentry.n_offset;
        else
	    entry->n_name = savename(symentry._n._n_name);
        entry->n_value = symentry.n_value;
        entry->n_scnum = symentry.n_scnum;
        entry->n_type = symentry.n_type;
        entry->n_sclass = symentry.n_sclass;
        entry->n_numaux = symentry.n_numaux;
        dprintf( (reloc) ? 5 : 9 ,
		"   entering symbol %-20s v=0x%x sn=%d ty=%d sc=%d na=%d ",
	            entry->n_name, entry->n_value, entry->n_scnum, 
                    entry->n_type,  entry->n_sclass, entry->n_numaux);

        switch (entry->n_scnum) {
            case N_DEBUG:   /* no need to enter debugging symbols now */
            case N_UNDEF:   /* no need to enter undefined externals */
	         dprintf( (reloc) ? 5 : 9 , "*\n");
                 break;
            case N_ABS:     /* absolute symbol - no need for reloc */
	         dprintf( (reloc) ? 5 : 9 , "\n");
                 enter(entry);
                 break;
            default:        /* some section number */
	         dprintf( (reloc) ? 5 : 9 , "\n");
	         if (reloc) entry->n_value += (long) textsegment;
	         enter(entry);
                 break;
        } /* switch */
    } /* for */
} /* readsymbols */



/***************************************************************/
/* processrelocation(nreloc, segment);                         */
/*                                                             */
/*     Processes the relocation information contained in the   */   
/* next nreloc bytes of the input file and makes the necessary */
/* adjustments to the memory beginning at the segment pointer. */
/***************************************************************/

static void processrelocation(nreloc, segment)
unsigned long nreloc;
char *segment;
{
#if 0
    struct relocation_info item;		NOT IN SYS5 */
#endif
    struct reloc item;                      /* SYS5 part */
    string name;
    symbol extsym;
    long offset;
    int i;


    dprintf(5,"processrelocation: nreloc = %d start @ 0x%x\n",
               nreloc,FTELL(ldptr));
    dprintf(5,"   sizeof item= %d   RELSZ = %d\n",
			sizeof(item), RELSZ);
    for (i = 0; i < nreloc; i++) {
/*	FREAD((char *) &item, sizeof item, 1, ldptr);	*/
	FREAD((char *) &item, RELSZ, 1, ldptr);	
        name = localtable[item.r_symndx].n_name;
	dprintf(5,"   looking for %-20s type=0x%x symndx=0x%x ",
			name,item.r_type,item.r_symndx);

            /* NEXT SECTION ONLY WORKS FOR EXTERN SYMBOLS - CAREFUL RECODE */
        extsym = lookup(name);
        if (extsym == NULL)
	    error("loadobj: undefined symbol %s\n", name);
	offset = extsym->n_value;
	dprintf(5," offset=0x%x\n",offset);
            /* ELSE (i.e. if local = static variable ) */
            /* offset = textsegment; */

	switch (RELOC_RTYPE(item)) {
            case R_POS:
                *((long *) &segment[item.r_vaddr]) += offset; break;
            case R_TOC:
                warning("R_TOC relocation not implemented yet");break;
            case R_BR:
                warning("R_BR relocation not implemented yet");break;
            case R_RELLONG:
                *((long *) &segment[item.r_vaddr]) += offset; break;
            case R_RELBYTE:
                segment[item.r_vaddr] += offset; break;
            default:
                error("loadobj: this type of reloc (0x%x) not yet supported",
                        RELOC_RTYPE(item));
	}
    }
}



/***************************************************************/
/* Symbol table functions:                                     */
/*      sym = lookup("str")  -- looks up string value in tree  */
/*      sym = enter(sym)     -- enters copy of sym in tree     */
/***************************************************************/

static symbol lookup(name)
string name;
{
    return (findentry(name, &symbase, (symbol) NULL));
}

static symbol enter(sym)
symbol sym;
{
    return (findentry(sym->n_name, &symbase, sym));
}



/***************************************************************/
/* Local routine that does the work for lookup and enter.      */
/* This is separated out because this level represents the     */
/* appropriate recursive formulation.                          */
/***************************************************************/

static symbol findentry(name, eptr, value)
string name;
symtree *eptr;
symbol value;
{
    int cmp;
    symtree entry;
    string savename;

    if ((entry = *eptr) == NULL) {      /* if at end a of branch */
	if (value == NULL) return NULL;     /* if looking -- nothing found */
	*eptr = entry = (symtree) getmem(sizeof *entry);  /* if entering */
	entry->stdata.n_name = scopy(value->n_name);    /* allocate and copy */
	entry->before = entry->after = NULL;
	cmp = 0;
    } else {
	cmp = strcmp(name, entry->stdata.n_name);       /* still checking */
    }
    if (cmp == 0) {                 /* found one or initializing */
	if (value != NULL) {   /* entering new symbol or overwriting old one */
	    savename = entry->stdata.n_name;
	    entry->stdata = *value;
	    entry->stdata.n_name = savename;
	}
	return (&entry->stdata);
    }
    return (findentry(name, (cmp<0) ? &entry->before : &entry->after, value));

}

/***************************************************************/
/* Local routine that saves a short (<= 8 char) symbol name in */
/* safe memory, as our information in 'nlist' items can only   */
/* have pointers to the name. Long names live in the string-   */
/* table which was already declared in safe memory in          */
/* readstrings()                                               */
/***************************************************************/

static string savename(shortname)
char shortname[8];
{
    int len;
    string cp;
    
    len = strlen(shortname);
    if (len>8)
        len=8;
    cp = getmem(len+1);
    strncpy(cp,shortname,len);
    *(cp+len) = '\0';
    return(cp);
}


#ifdef TESTBED_OLD

#include <getparam.h>

string defv[] = {
    "func=(n==0?1:n*f(n-1))\n   Function to test, or to call",
    "low=1\n                    Low loop index for test run",
    "high=10\n                  High loop index for test run",
    "object=\n                  Name of object file to load",
    "type=\n                    Function type (r,i,v)",
     NULL,
};

double sin(), cos(), sqrt();
local int nmain = 0;

main(argc, argv)
int argc;
string argv[];
{
    string fdef, object, type;
    int    low, high, i;
    stream tmpfile;
    double x;
    proc   fn;
    iproc  fni;
    rproc  fnr;

    initparam(argv, defv);

    fdef = getparam("func");
    low = getiparam("low");
    high = getiparam("high");
    object = getparam("object");
    type = getparam("type");

    mysymbols(getparam("argv0"));

    if (*object) {
        loadobj(object);         /* load it */
        if (*fdef) {            /* if a function name supplied ... */
            switch (*type) {    /* determine it's type */
              case 'r':
                    fnr = (rproc) findfn(fdef);  /* find name */
                    if (fnr==NULL) error("Error finding %s\n",fdef);
                    for (i=low; i<=high; i++)
                        (*fnr)((double)i);          /* and call it */
                    break;
              case 'i':
                    fni = (iproc) findfn(fdef);
                    if (fni==NULL) error("Error finding %s\n",fdef);
                    for (i=low; i<=high; i++)
                        (*fni)(i);
                    break;
              default:
                    fn = (proc) findfn(fdef);
                    if (fn==NULL) error("Error finding %s\n",fdef);
                    for (i=low; i<=high; i++)
                        (*fn)();
                    break;
            } /* switch */
        } /* if(func) */
        exit(0);
    } /* if(object) */

    if (*fdef != '0') {
    tmpfile = fopen("ld-tmp.c", "w");               /****** TMPFILE ******/
    fprintf(tmpfile, "#include <stdio.h>\n");
    fprintf(tmpfile, "double sin(), cos(), tryext();\n");
    fprintf(tmpfile, "extern int nmain;\n");
    fprintf(tmpfile, "static int nlocal=0;\n");
    fprintf(tmpfile, "int f(n)\n");
    fprintf(tmpfile, "int n;\n");
    fprintf(tmpfile, "{\n");
    fprintf(tmpfile, "    dprintf(1,\" calling f(n=%%d)\\n\",n);\n");
    fprintf(tmpfile, "    nmain += n;\n");
    fprintf(tmpfile, "    tryext((double)n);\n");	
    fprintf(tmpfile, "    return (%s);\n", fdef);
    fprintf(tmpfile, "}\n");
    fclose(tmpfile);                                /**********************/
    if (system("cc -c ld-tmp.c") != 0)             
	error("function %s does not parse\n", fdef);
    } else
	printf("Attempt to load local file ld-tmp.o\n");
    loadobj("ld-tmp.o");
#if 0
    if (system("rm -f ld-tmp.c ld-tmp.o") != 0)
	error("cannot rm ld-tmp.*\n");
#endif
    fn = (iproc) findfn("f");	        /* SYS5 has no leading underscore */
    if (fn == NULL)
	error("function not correctly defined\n");
    for (i = low; i <= high; i++) {
	printf("f(%d) = %d,  nmain=%d\n", i, (*fn)(i),nmain);
        x++;
    }
    exit(0);
}

double tryext(x)
double x;
{
  printf("Tryext (%g) called\n",x);
  return (x+1.0);
}

_math_loader()		/* ensures loading of some math functions */
{			/* this routine itself is never called    */
    sin(1.0);
    cos(1.0);
    sqrt(1.0);
}
#endif
