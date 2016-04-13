#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif



int check_for_output ( char *search_string );
void send_command(char * sent_command);
double re_estimate_parameters (double likelihood);
double build_starting_tree (double likelihood);
int nexusparser(FILE *nexusfile);
int comment(FILE *file);
int setup_pipe (void);
int close_pipe (void);
void unroottree(char * tree);
int test_reverse_constraints(char * translated_tree);
double build_starting_constraint_tree (double likelihood, char *constraint, int constraint_num); 
double re_estimate_constraint_parameters (double likelihood, int constraint_num);




FILE *paup_pipe, *logfile, *treefile, *outputtree;

int bytes_read, pid,  file_length=0, checks=0, checkcount=0;
int nbytes = 100, numtranslatedtaxa = 0, num_taxa = 1000, translated = FALSE, alignment_length =0, number_taxa =0;
char *my_string, pid_string[1000], c, *newtree, *notranslate_newtree;
char logfilename[1000], sys[1000], result[100000], output[100000], command[10000], *treefilename, ***names;
double likelihood =0, new_likelihood=0, meanRand_likelihood=0;




int main (int argc, char *argv[])
{
  char * token, string[1000000];
  int iteration=1, error = FALSE, i, j;
  
  newtree = malloc(100000*sizeof(char));
  notranslate_newtree = malloc(100000*sizeof(char));
  newtree[0] = '\0'; notranslate_newtree[0] = '\0';

  pid=getpid();
  sys[0] = '\0';
  pid_string[0]='\0';
  sprintf(pid_string, "%d", pid);
  command[0] = '\0';
  result[0]='\0';

  if(argc < 2)
  {
    printf("\n\nMachete: Likelihood reverse constraint analysis using PAUP\n\n Usage: \"machete <nexus file>\"\n\nWhere: <nexus file> is a nexus formatted alignment file of DNA sequences\n\n" );
    exit(1);
  }


  /* initiate PAUP pipe */
  if(setup_pipe () == TRUE)
    {

    printf("\nExecuting file %s\n", argv[1]);
    sprintf(command, "execute %s;\n", argv[1]);
    send_command (command);
    checkcount =0; while(!check_for_output("Data matrix has")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    token = strtok(output, " ");
    for(i=0; i<3; i++) token = strtok(NULL, " ");
    number_taxa = atoi(token);
    for(i=0; i<2; i++) token = strtok(NULL, " ");
    alignment_length = atoi(token);
    printf("\tNumber of taxa = %d Number of Characters = %d\n\n", number_taxa, alignment_length);


    likelihood = build_starting_tree (likelihood);

    likelihood = re_estimate_parameters (likelihood);

    /* now read tree from PAUP file, capturing the name definitions */
    treefilename = malloc(1000*sizeof(char));
    treefilename[0] = '\0';
    checkcount =0; while(!check_for_output("saved to file")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    sprintf(treefilename, "paup_%d.tre", pid);
    treefile = fopen(treefilename, "r");
    nexusparser(treefile);
    fclose(treefile);
    if(!translated) strcpy(notranslate_newtree, newtree);
    unroottree(newtree);
    unroottree(notranslate_newtree);


    /* carry out randomisation test to identify mean of random trees */
    /*printf("Carrying out random trees analysis\n");
    send_command("randtrees;\n");
    checkcount =0; while(!check_for_output("mean=")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}

    token = strtok(output, " =");
    token = strtok(NULL, " =");

    meanRand_likelihood=atof(token);
    printf( "\n\tMean -ln L of random trees: %g\n", meanRand_likelihood );
*/
/*
    printf("%s\n", newtree);
    
    printf("%s\n", notranslate_newtree);

     if(translated)
      {
      for(i=0; i<numtranslatedtaxa; i++)
        {
        for(j=0; j<2; j++)
          printf("%s\n", names[i][j]);
        }
      }

*/

    /* we now have the best tree in memory anlong with a translation table */

    /* Next identify each split and for each, define a constraint, and fine the likelihood of the best tree that does not contain that split */
      printf("\nExtract all constraints defined by best tree and carry out reverse constraint analysis for each\n");

     test_reverse_constraints(newtree);

     

     sprintf(command, "%s.labelledtree.tre", argv[1]);
     outputtree = fopen(command, "w");
     fprintf(outputtree, "%s [%g]\n", newtree, likelihood);
     fclose (outputtree);

     printf("\n\n\nresulting weighted tree:\n%s [%g]\n\nTree with labels saved to file %s\n\n", newtree, likelihood, command);

   /* Clean up and exit */
    free(treefilename);

    /*unroottree(string); */

    error = close_pipe();
    
    }
  else
    {
    error = TRUE;
    }
  free(newtree);

  if(translated)
    {
    for(i=0; i<num_taxa; i++)
      {
      for(j=0; j<2; j++)
        free(names[i][j]);
      free(names[i]);
      }
    free(names);
    }


return(error);
}







int check_for_output ( char *search_string )
  {
  int found=FALSE, i=0;


 /* printf("searching for %s\n", search_string); */
  
  /* look for a specific output from the command */
  output[0] ='\0';

  sleep(1);
  logfile=fopen(logfilename, "r"); 
   
  if(logfile == '\0') 
    {
      printf("Problem opening logfile\n!");
    }
  else
    {   
    i=0;
    while(!feof(logfile) && !found) {
          output[i]= getc(logfile);
          i++;
          if(output[i-1] == '\n' || output[i-1] == '\r')
          {
          output[i] = '\0';
          if(strstr(output, search_string ) != '\0')
              found = TRUE;
          i=0;
          }
      }
    }
  fclose(logfile);

  /*** End looking for output */
  return(found);
  }



void send_command(char * sent_command)
  {
  char com[1000];

  com[0] = '\0';

  /* delete logfile, so that outputs from previous commands do not get mixed up with the outputs of the command about to be run */
  remove(logfilename);
  strcpy(com, "touch ");
  strcat(com, logfilename);
  system(com);

 /* printf ("Command sent: %s\n", sent_command); */
  fprintf (paup_pipe, "log file=%s start=yes replace=yes flushLog=yes; \n", logfilename);
  fflush(paup_pipe);

  fprintf (paup_pipe, "%s\n", sent_command);
  fflush(paup_pipe);
  fprintf (paup_pipe, "log stop=yes;\n");
  fflush(paup_pipe);
  sleep(2);
  }


double re_estimate_parameters (double likelihood)
  {
  double new_likelihood = 0;
  int i=0, iteration=1;
  char * token, *comm;

  comm = malloc(10000*sizeof(char));
  comm[0] = '\0';
   while(round(new_likelihood) < round(likelihood))
      {
      printf("Iteration %d: Estimating parameters given tree and recalculating likelihood", iteration);

      send_command ("lset Rmatrix=estimate Basefreq=Estimate rates=gamma shape=estimate; lscores;\n");
     
      checkcount=0; while(!check_for_output("-ln L")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
      token = strtok(output, " ");
      for(i=0; i<2; i++) token = strtok(NULL, " ");
      new_likelihood=atof(token);
      printf( "\n\tBest -ln L: %g\n", new_likelihood );

      send_command ("lset Rmatrix=previous Basefreq=previous rates=gamma shape=previous; hs;");

      if(round(new_likelihood) < round(likelihood))
        {
        printf("Iteration %d: Better likelihood found, Saving new estimated parameters. Beginning tree search", iteration);
        likelihood = new_likelihood;
        checkcount=0; while(!check_for_output("Score of best tree(s) found")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
        token = strtok(output, " ");
        for(i=0; i<6; i++) token = strtok(NULL, " ");
        new_likelihood=atof(token);
        printf( "\n\tBest -ln L: %g\n", new_likelihood );
        }
      iteration++;
      }
    printf("\tDifference in -ln L less than 1 unit, stopping search\n");
    sprintf(comm, "savetrees file=paup_%d.tre brLens=yes;\n", pid);
    send_command(comm);
    printf("Best tree found after %d iterations with -ln L %g saved to file paup_%d.tre;\n\n", iteration-1,  new_likelihood, pid); 

    free(comm);
    return(new_likelihood);
  }


  double build_starting_tree (double likelihood)
    {
    char * token;
    int i;

    printf("Generating starting ML tree\n");
    send_command ("set criterion=like; lset NST=6; hs;");
 
    checkcount=0; while(!check_for_output("Score of best tree(s) found")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    token = strtok(output, " ");
    for(i=0; i<6; i++) token = strtok(NULL, " ");
    likelihood=atof(token);
    printf( "\tBest -ln L: %g\n\n", likelihood );

    return(likelihood);
    }

  


 int setup_pipe (void)
 { 
  int success = TRUE;

  /* Open the pipe */
  paup_pipe = popen ("paup -n > Paup_output.txt 2> /dev/null", "w");

  /* Check that pipes are non-null, therefore open */
  if ((!paup_pipe))
    {
      fprintf (stderr,
               "Paup pipe failed.\n");
      success = FALSE;
    }
  else
    {  

    /* Stop Paup from printing output unnecessarily to the screen */
    fprintf (paup_pipe, "set increase=auto warnReset=no warnTree=no warnTSave=no warnBlkName=no warnRoot=no  warnRedef=no errorStop=no errorBeep=no queryBeep=no flushLog=yes;\n");
 
    
    /* Create logfilename for Paup */
    logfilename[0]='\0';
    strcpy(logfilename, "paup_");
    strcat(logfilename, pid_string);
    strcat(logfilename, "_log.txt");
   /* printf("%s\n", logfilename);*/

    /* create the log file - so it doesn;t cause an error when Machete goes looking for it */
    strcpy(sys, "touch ");
    strcat(sys, logfilename);
    system(sys);

    }

  return(success);
 }


int close_pipe ()
  {
  int success = TRUE;
   /* Quit Paup */
    fprintf (paup_pipe, "quit;\n");
    fflush(paup_pipe);

    /* Close grep_pipe, cehcking for errors */
    if (pclose (paup_pipe) != 0)
      {
        fprintf (stderr,
                 "Could not run 'paup', or other error.\n");
        success = FALSE;
      }

  return(success);
  }




int nexusparser(FILE *nexusfile)
  {
  char c, begin[6] = {'b', 'e', 'g', 'i', 'n', '\0'}, translate[10] = {'t','r','a','n','s','l','a','t','e','\0'}, tree[5] = {'t','r','e','e','\0'};
  int error = FALSE, i, j, k, l, x,  num_trees = 0, found = FALSE;
  char *string, single[1000];
  

  string = malloc(400000*sizeof(char));
  string[0] = '\0';
  while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));

  while(!feof(nexusfile) && !error)
    {
    while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
    switch(tolower(c))
      {
      case '[':
        
        comment(nexusfile);
        break;
      
      case 'b':   /* this is the beginning of a block **/
        for(i=1; i<5; i++)
          if((c = tolower(getc(nexusfile))) != begin[i]) error = TRUE;
        if(!error)
          {
          while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
          i = 0;
          do {  /* find the type of block it is **/
            string[i] = tolower(c);
            i++;
            c = getc(nexusfile);
            }while(c != ';' && c != ' ');
          string[i] = '\0';
          if(strcmp(string, "trees") == 0)
            {
            while(strcmp(string, "end;") != 0 && strcmp(string, "endblock;") != 0 && !feof(nexusfile) && !error)
              {
              /*** read in the trees ***/
              while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
              while(c == '[')
                {
                comment(nexusfile);
                while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
                }
              strcpy(string, "");
              i=1;
              string[0] = tolower(c);
              while((c = getc(nexusfile)) != ' ' && c != '\n' && c != '\r' && !feof(nexusfile))
                {
                string[i] = tolower(c);
                i++;
                }
              string[i] = '\0';
              if(strcmp(string, tree) == 0)
                {
                /** read in trees **/
                /** read in the name of the tree **/
                while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
                i=1;
                string[0] = c;
                while(((c = getc(nexusfile)) != ' ' && c != '\n' && c != '\r') && !feof(nexusfile))
                  {
                  string[i] = c;
                  i++;
                  }
                string[i] = '\0';
                /*strcpy(tree_names[Total_fund_trees], string); */ /*  not needed for this implementation */
                
                /*** look out for any weights **/
                while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '=' || c == '\n' || c == '\r') && !feof(nexusfile));
                while(c == '[')
                  {
                  if((c = getc(nexusfile)) == '&')
                    {
                    if((c = tolower(getc(nexusfile))) == 'w')
                      {
                      while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
                      i=1;
                      string[0] = c;
                      while((c = getc(nexusfile)) != ']' &&!feof(nexusfile))
                        {
                        string[i] = c;
                        if(string[i] != ' ') i++;
                        }
                      string[i] = '\0';
                     /* tree_weights[Total_fund_trees] = tofloat(string); */ /* Not needed for this implementation */
                      }
                    else
                      while((c = getc(nexusfile)) != ']' &&!feof(nexusfile));
                    }
                  else
                    while((c = getc(nexusfile)) != ']' &&!feof(nexusfile));
                  while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r') && !feof(nexusfile));
                  }
                
                
                /*** read in the tree **/
                if(!translated)
                  {
                  i=0;
                  strcpy(string, "");
                  while(c != ';')
                    {
                    if(string[i] != ' ' && string[i] != '\t' && string[i] != '\n' && string[i] !='\r') string[i] = c;
                    i++;    
                    c = getc(nexusfile);
                    }
                  string[i] = ';';
                  string[i+1] = '\0';
                 /* input_fund_tree(string, Total_fund_trees); */ /* Not needed for this implementation */
                  }
                else
                  {
                  /*** If there was a translation table then we need to put the names into the tree **/
                  i=0; x=0;
                  while(c != ';' && !feof(nexusfile) && !error)
                    {
                    switch(c)
                      {
                      case '(':
                      case ',':
                        newtree[i] = c;
                        notranslate_newtree[x] = c;
                        i++; x++;
                        c = getc(nexusfile);
                        break;
                      case ')':
                      case ':':
                        do
                          {
                          newtree[i] = c;
                          notranslate_newtree[x] = c;
                          i++; x++;
                          c = getc(nexusfile);
                          }while(c != '(' && c != ')' && c != ',' && c != ';');
                        break;
                      default:
                        j=0;
                        while(c != '(' && c != ')' && c != ',' && c != ';' && c != ':' && !feof(nexusfile))
                          {
                          single[j] = c;
                          notranslate_newtree[x] = c;
                          j++; x++;
                          c = getc(nexusfile);
                          }
                        single[j] = '\0';
                        found = -1;
                        for(j=0; j<numtranslatedtaxa; j++)
                          {
                          if(strcmp(single, names[j][0])== 0)
                            found = j;
                          }
                        if(found == -1)
                          {
                          printf("Error: taxa %s is not in the translation table\n", single);
                          error = TRUE;
                          }
                        else
                          {
                          j=0;
                          while(names[found][1][j] != '\0')
                            {
                            newtree[i] = names[found][1][j];
                            j++; i++;
                            }
                          }
                      }   
                    }
                  newtree[i] = ';';
                  notranslate_newtree[x] = ';';
                  newtree[i+1] = '\0'; notranslate_newtree[x+1] = '\0';
                 /* input_fund_tree(newtree, Total_fund_trees); */  /* Not needed for this implementation */
                  }
                }
              if(strcmp(string, translate) == 0)
                {
                translated = TRUE;
                /** handle translation **/
                names = malloc(num_taxa*sizeof(char**));
                for(i=0; i<num_taxa; i++)
                  {
                  names[i] = malloc(2*sizeof(char*));
                  for(j=0; j<2; j++)
                    {
                    names[i][j] = malloc(1000*sizeof(char));
                    names[i][j][0] = '\0';
                    }
                  }
                i=0; j=0; c = getc(nexusfile);
                while(c != ';' && !feof(nexusfile))
                  {
                  while(c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == ',') c = getc(nexusfile);
                  if(j==2)
                    {
                    j = 0;
                    i++;
                    if(i+1 >= num_taxa)
                      {
                      names = realloc(names,(num_taxa+num_taxa)*sizeof(char**));
                      for(k=num_taxa; k<num_taxa+num_taxa; k++)
                        {
                        names[k] = malloc(2*sizeof(char*));
                        for(l=0; l<2; l++)
                          {
                          names[k][l] = malloc(1000*sizeof(char));
                          names[k][l][0] = '\0';
                          }
                        }
                      
                      num_taxa += num_taxa;
                      }
                    }
                  k=0;
                  while(c != ' ' && c != '\n' && c != '\r' && c != '\t' && c != ';' && c != ',' && !feof(nexusfile))
                    {
                    names[i][j][k] = c;
                    k++;
                    c = getc(nexusfile);
                    }
                  names[i][j][k] = '\0';
                  j++;
                  }
                numtranslatedtaxa = i+1;
                }
              }
            }
          else
            {
            /*** skip to the end of the block -- it must contain sequences or something **/
            strcpy(string, "");
            while(strcmp(string, "end;") != 0 && strcmp(string, "endblock;") != 0 && !feof(nexusfile) && !error)
              {
              i=0;
              while((c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == ';') && !feof(nexusfile)) c = getc(nexusfile);
              while(c != ';' && c != ' ' && c != '\n' && c != '\r' && c != '\t' && !feof(nexusfile))
                {
                string[i] = tolower(c);
                if(i < 9997)i++;
                c = tolower(getc(nexusfile));
                }
              string[i] = ';';
              string[i+1] = '\0';
              }
            } 
          }
        else
          {
          /*** what to do if it doesn't say begin? ***/
          error = TRUE;
          }
        break;
        
      default:
        break;


      }
    }
  
/*  fclose(nexusfile); */
  free(string);
  return(error);
  }


int comment(FILE *file)
  {
  char c;
  c = getc(file);
  while(c != ']' && !feof(file))
    {
    if(c == '[') comment(file);
    c = getc(file);
    }
  if(!feof(file)) return(FALSE);
  else return(TRUE);
  }


void unroottree(char * tree)
    {
    int i=0, j=0, k=0, l=0, m=0, basecount = 0, parentheses=0;
    int foundopen = FALSE, foundclose = FALSE;
  float del_nodelen = 0;
  char length[100], restof[400000];
  
  restof[0] = '\0';
  length[0] = '\0';
  /* scan through the tree counting the number of taxa/nodes at the base (for it to be unrooted there should be at least three) */
    while(tree[i] != ';')
        {
        switch(tree[i])
            {
            case '(':
                    parentheses++;
                    i++;
                    break;
            case ')':
                    parentheses--;
                    i++;
                    break;
            case ',':
                    if(parentheses == 1)
                        {
                        basecount++;
                        }
                    i++;
                    break;
            default:
                    i++;
                    break;
            }
        }
        
    if(basecount <2)  /* if the base of the tree is rooted */
        {
        i=0;
        parentheses = 0;
        while(tree[i] != ';')  /* delete the two parentheses to make the tree unrooted */
            {
            switch(tree[i])
                {
                case '(':
                        parentheses++;
                        if(parentheses == 2 && !foundopen)
                            {
                            tree[i] = '^';
                            foundopen = TRUE;
                            }
                        i++;
                        break;
                case ')':
                        if(parentheses == 2 && !foundclose)
                            {
                            tree[i] = '^';
              i++;
              while(tree[i] != ')' && tree[i] != '(' && tree[i] != ',' && tree[i] != ';' && tree[i] != ':')
                {
                tree[i] = '^';
                i++;
                }
                            if(tree[i] == ':')
                                {
                k=0;
                length[0] = '\0';
                                while(tree[i] != ')' && tree[i] != '(' && tree[i] != ',' && tree[i] != ';')
                                    {
                  if(tree[i] != ':')
                    {
                    length[k] = tree[i];
                    k++;
                    }
                  tree[i] = '^';
                                    i++; 
                                    }
                length[k] = '\0';
                                }
              if(length[0] != '\0') /* we have a branch length on the internal branch, so we need to add it to the branch length of the component that is to the direct right of this closing parenthesis */
                {
                del_nodelen = atof(length);
                k=i+1; /* This should be whatever is after the ',' which should be the next compoonent */
                if(tree[k] == '(') /* we need to find the end of this clade and add the value there */
                  {
                  l=1; k++;
                  while((l != 0 || tree[k-1] != ')') && tree[k] != ';' )   /* CHANGED RECENTLY FROM while(l != 0 && tree[k-1] != ')' && tree[k] != ';' ) */
                    {
                    switch(tree[k])
                      {
                      case '(':
                        l++;
                        k++;
                        break;
                      case ')':
                        l--;
                        k++;
                        break;
                      default:
                        k++;
                        break;
                      }
                    }
                  k--; /* k now points to the closing bracket */
                  /* read in the length attached to this partenthsis */
                  k++;
                  while(tree[k] != ')' && tree[k] != '(' && tree[k] != ',' && tree[k] != ';' && tree[k] != ':') k++; /* k now points to the ":" at the start of the length (if there is a length) */
                  }
                else
                  {
                  while(tree[k] != ')' && tree[k] != '(' && tree[k] != ',' && tree[k] != ';' && tree[k] != ':') k++; /* k now points to the ":" at the start of the length (if there is a length) */
                  }
                if(tree[k] == ':') /* there is length attached to this */
                  {
                  m=k+1;
                  length[0] = '\0'; l=0;
                  while(tree[m] != ')' && tree[m] != '(' && tree[m] != ',' && tree[m] != ';')
                    {
                    length[l] = tree[m];
                    l++; m++;
                    }
                  length[l] = '\0';
                  del_nodelen += atof(length);
                  
                  }
                else
                  m=k;
                /* now add del_nodelen to this point in the tree */
                 l=0;
                while(tree[m] != ';' && tree[m] != '\0')
                  {
                  restof[l] = tree[m];
                  m++; l++;
                  }
                restof[l] = ';';
                restof[l+1] = '\0';
                if(tree[k] == ':')
                  tree[k] = '\0';
                else
                  {
                  tree[k] = '\0';
                  }
                length[0] = '\0';
                sprintf(length, ":%f", del_nodelen);
                strcat(tree, length);
                strcat(tree, restof);
                }
                            foundclose = TRUE;
                            }
                        i++;
                        parentheses--;
                        break;
                default:
                        i++;
                        break;
                }
            }
                        
        /* scan through the string shifting up the characters to take into account those parentheses that have been deleted */
        i=0; j=0;
        while(tree[j] != ';')
            {
            if(tree[j] == '^')
                {
                while(tree[j] == '^')
                    j++;
                if(i!= j)tree[i] = tree[j];
                i++; j++;
                }
            else
                {
                if(i!=j)tree[i] = tree[j];
                i++;j++;
                }
            }
        tree[i] = tree[j];
        tree[i+1] = '\0';
        
        }
    }


int test_reverse_constraints(char * translated_tree)
  {
    int i=0, k=0, l=0, j=0, x=0, y=0, z=0, q=0, splitnumber=0, constraint_num=0;
    char taxaname[1000], constraint[1000000], weight[10000], tree[1000000], tmptree[1000000];

    taxaname[0] = '\0';
    constraint[0] = '\0'; j=0;
    weight[0] = '\0';
    tmptree[0] = '\0';
    tree[0]='\0';

    

    strcpy(tree, translated_tree);

    while(tree[i] != ';')
    {
      if( tree[i] == '(' && i !=0)
        {
          /* identify all the taxa IDs in this split */
          k=i; j=0;
          constraint[j] = '('; j++;
           l=1; k++;
          while((l != 0 || tree[k-1] != ')') && tree[k] != ';' )   /* CHANGED RECENTLY FROM while(l != 0 && tree[k-1] != ')' && tree[k] != ';' ) */
            {
            switch(tree[k])
              {
              case '(':
                l++;
                k++;
                break;
              case ')':
                l--;
                k++;
                break;
              case ':':
                while(tree[k] != '(' && tree[k] != ')' && tree[k] != ',' && tree[k] != ';' ) k++;
                break;
              case ',':
                k++;
                break;
              default:
                while(tree[k] != ',' && tree[k] != '(' && tree[k] != ')' && tree[k] != ':' )
                    {
                    constraint[j] = tree[k];
                    k++; j++;
                    }
                constraint[j] = ',';
                j++;
                break;
              }
            }
          constraint[j-1] = ')'; /* overwrites the extra comma at the end */
          constraint[j] = '\0';
          

          /* Calculate the likeliood of the best tree that does not contain this constraint */


          new_likelihood = build_starting_constraint_tree (likelihood, constraint, constraint_num);

          new_likelihood = re_estimate_constraint_parameters (new_likelihood, constraint_num);


          printf("\tBest reverse Constraint -ln L = %g\tdifference = %g\n", new_likelihood, new_likelihood-likelihood);
         /* sprintf(weight, "%g/%g", new_likelihood- likelihood, (new_likelihood- likelihood)/(meanRand_likelihood-likelihood)); */
          sprintf(weight, "%g", new_likelihood- likelihood);


          /* find the end of this split in the named tree string to append the weight as a label */

          /* find the beginning of this split in the translated tree */
          x=0; y=-1;
          while(translated_tree[x] != '(' || y!=constraint_num)
            {
            x++;
            if(translated_tree[x] == '(') y++;
            }
          /* x now points to the beginning of the current split */
          /* find the end of the split */
          l=0;
         while((l != 0 || translated_tree[x-1] != ')') && translated_tree[x] != ';' )   /* CHANGED RECENTLY FROM while(l != 0 && tree[k-1] != ')' && tree[k] != ';' ) */
            {
            switch(translated_tree[x])
              {
              case '(':
                l++;
                x++;
                break;
              case ')':
                l--;
                x++;
                break;
              default:
                x++;
                break;
              }
            }
          /* x now points to the position AFTER the closing bracket for this component */
            q=0;
            for(z=x; z<(strlen(translated_tree)); z++)
              {
              tmptree[q] = translated_tree[z];
              q++;
              }
            tmptree[q] = '\0';

          /* add weight to end of component */
            translated_tree[x] = '\0';
            strcat(translated_tree, weight);
            strcat(translated_tree, tmptree);



          constraint_num++;
        }
     i++;   
    }

    return(1);
  }




double build_starting_constraint_tree (double likelihood, char *constraint, int constraint_num)
    {
    char * token, *comm;
    int i;


    comm = malloc(10000*sizeof(char));
    comm[0] = '\0';
    printf("\t\n====================\nSetting constraint number %d: %s\n====================\n", constraint_num, constraint);
    sprintf(comm, "constraints constraint_%d = (%s);", constraint_num, constraint);
    send_command(comm );

    printf("\tGenerating starting ML tree using reverse constraint tree\n");
    sprintf(comm, "hs constraints=constraint_%d enforce=yes converse=yes;\n", constraint_num);
    send_command (comm);
 
    checkcount=0; while(!check_for_output("Score of best tree(s) found")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    token = strtok(output, " ");
    for(i=0; i<6; i++) token = strtok(NULL, " ");
    likelihood=atof(token);
    printf( "\t\tBest -ln L: %g\n\n", likelihood );

    free(comm);
    return(likelihood);
    }




double re_estimate_constraint_parameters (double likelihood, int constraint_num)
  {
  double new_likelihood = 0;
  int i=0, iteration=1;
  char * token, *comm;

  comm = malloc(10000*sizeof(char));
  comm[0] = '\0';
   while(round(new_likelihood) < round(likelihood))
      {
      printf("\tIteration %d: Estimating parameters given tree and recalculating likelihood", iteration);

      send_command ("lset Rmatrix=estimate Basefreq=Estimate rates=gamma shape=estimate; lscores;\n");
     
      checkcount=0; while(!check_for_output("-ln L")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
      token = strtok(output, " ");
      for(i=0; i<2; i++) token = strtok(NULL, " ");
      new_likelihood=atof(token);
      printf( "\n\t\tBest -ln L: %g\n", new_likelihood );


      send_command ("lset Rmatrix=previous Basefreq=previous rates=gamma shape=previous;\n");

      if(round(new_likelihood) < round(likelihood))
        {
        printf("\tIteration %d: Better likelihood found, Saving new estimated parameters. Beginning tree search using reverse constraint", iteration);

        likelihood = new_likelihood;
        

        sprintf(comm, "hs constraints=constraint_%d enforce=yes converse=yes;\n", constraint_num);
        send_command (comm);
 

        checkcount=0; while(!check_for_output("Score of best tree(s) found")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
        token = strtok(output, " ");
        for(i=0; i<6; i++) token = strtok(NULL, " ");
        new_likelihood=atof(token);
        printf( "\n\t\tBest -ln L: %g\n", new_likelihood );
        }
      iteration++;
      }
    printf("\t\tDifference in -ln L less than 1 unit, stopping search\n");

    free(comm);
    return(new_likelihood);
  }
