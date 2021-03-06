#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
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


time_t interval1, interval2, interval3;
double telapsed, diff;
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
void print_time(double diff);
void read_sitelike_file(int component);
 void do_bootstrap_analysis (int reps);






FILE *paup_pipe, *logfile, *treefile, *outputtree, *nexusfile, *guidetreefile;

int bytes_read, pid,  file_length=0, checks=0, checkcount=0, total_constraints=0, response=0, datatype=0, constart=0, conend=-1, bootreps=0; /*dataypes 0=DNA, 1=Protein */
int nbytes = 100, numtranslatedtaxa = 0, num_taxa = 1000, translated = FALSE, alignment_length =0, excludedchars=0, number_taxa =0, print_command=FALSE, gtflag=FALSE, num_trees=0, buildflag=FALSE, listcons=FALSE, norevcon=FALSE;
char *my_string, pid_string[1000], c, *newtree, *notranslate_newtree, *besttree, **resultingtrees = NULL, *guidetree = NULL,  *infile=NULL, *gtfilename=NULL ;
char logfilename[1000], sys[1000], result[100000], output[100000], command[10000], *treefilename, ***names;
double likelihood =0, new_likelihood=0, meanRand_likelihood=0, **site_likelihoods = NULL;




int main (int argc, char *argv[])
{
  char * token, string[1000000], *tok, *str;
  int iteration=1, error = FALSE, i, j, found=FALSE, fflag=0, deltmp=TRUE;
  FILE *file = NULL;
  
  newtree = malloc(100000*sizeof(char));
  besttree = malloc(100000*sizeof(char));
  guidetree = malloc(100000*sizeof(char));
  str = malloc(100000*sizeof(char));
  notranslate_newtree = malloc(100000*sizeof(char));
  newtree[0] = '\0'; besttree[0] = '\0'; notranslate_newtree[0] = '\0';

  pid=getpid();
  sys[0] = '\0';
  pid_string[0]='\0';
  sprintf(pid_string, "%d", pid);
  command[0] = '\0';
  result[0]='\0';

  if(argc < 2)
  {
    printf("\n\nMachete: Likelihood reverse constraint analysis using PAUP\n\n Usage: \"machete -f <nexus file> -[cthln] [-s INTEGER] [-e INTEGER] [-r INTEGER]\"\n\n\tWhere: <nexus file> is a nexus formatted alignment file of DNA sequences\n\t-c commands sent to Paup to be also printed to standard error\n\t-t preserves temporary files\n\t-h prints this message\n\t-b force build optimum tree (when a tree has been provided in the nexus file)\n\t-s <constraint number> specifies the constraint to start at\n\t-e <constraint number> specific the constraint to end at\n\t-l list constraints (and do not carry out reverse constraints analysis)\n\t-n tells machete NOT to carry out the reverse constraint analysis (Just build the best tree)\n\t-r tells machete how many boostrap replicates to carry out (by default = 0)\n\n" );
    printf("\tCreated by Chris Creevey e:chris.creevey@gmail.com t:@hairy_llama\n\n\n");
    exit(1);
  }
  while ((c = getopt(argc, argv, "f:chtbs:e:lnr:")) != -1)
    {   
      switch (c) 
      {
      case 'c':
        print_command = TRUE;
        break;
      case 'f':
        fflag = 1;
        infile = optarg;
        break;
      case 'h':
        printf("\n\nMachete: Likelihood reverse constraint analysis using PAUP\n\n Usage: \"machete -f <nexus file> -[cthln] [-s INTEGER] [-e INTEGER] [-r INTEGER]\"\n\n\tWhere: <nexus file> is a nexus formatted alignment file of DNA sequences\n\t-c commands sent to Paup to be also printed to standard error\n\t-t preserves temporary files\n\t-h prints this message\n\t-b force build optimum tree (when a tree has been provided in the nexus file)\n\t-s <constraint number> specifies the constraint to start at\n\t-e <constraint number> specific the constraint to end at\n\t-l list constraints (and do not carry out reverse constraints analysis)\n\t-n tells machete NOT to carry out the reverse constraint analysis (Just build the best tree)\n\t-r tells machete how many boostrap replicates to carry out (by default = 0)\n\n" );
        printf("\tCreated by Chris Creevey e:chris.creevey@gmail.com t:@hairy_llama\n\n\n");

        exit(1);
        break;
      case 't':
        deltmp=FALSE;
        break;
      case 'b':
        buildflag=TRUE;
        break;
      case 's':
        constart=atoi(optarg);
        break;
      case 'e':
        conend=atoi(optarg);
        break;
      case 'l':
        listcons=TRUE;
        break;
     case 'n':
        norevcon=TRUE;
        break;
     case 'r':
        bootreps=atoi(optarg);
        if(bootreps == 0) 
          {
            printf ("ERROR: number of bootstrap replicates to be performed must be greater than 0\n");
            exit(0);
          }
        break;
      }
    }

  if(fflag == 0) 
    { /* -f was mandatory */
    fprintf(stderr, "%s: missing -f option\n", argv[0]);
    exit(1);
    }


 /* Test input file exists and read in the nexus file so we can determine the datatype */
if((nexusfile = fopen(infile, "r")) == '\0')   /* check to see if the file is there */
    {                          /* Open the fundamental tree file */
    fprintf(stderr, "Error: Cannot open nexus file %s\n", infile);
    exit(1);
    }
  else
    {
    while(!feof(nexusfile) && found==FALSE) 
      {
      str[0]='\0';
      fscanf(nexusfile,"%s", str);
      if(strstr(str, "datatype") || strstr(str, "DATATYPE"))
        { 
        found=TRUE; 
        if((strstr(str, "protein")) || (strstr(str, "PROTEIN")) ) datatype=1;
        else
          { 
            if((strstr(str, "dna")) || (strstr(str, "DNA")) || (strstr(str, "nucleotide")) || (strstr(str, "NUCLEOTIDE")) || (strstr(str, "rna")) || (strstr(str, "RNA"))) datatype=0;
            else 
              {
                printf("Error: Datatype in nexus file not recognised. It muct be one of \"dna\", \"rna\", \"nucleotide\" or \"protein\"\n");
                fclose(nexusfile);
                exit(0);
              }
          }
        }
      }
    if(found == FALSE)
      { 
      fprintf(stderr, "Error: Datatype not found in nexus file. It muct be one of \"dna\", \"rna\", \"nucleotide\" or \"protein\"\n");
      fclose(nexusfile);
      exit(0);
      }
    else
      {
      if(datatype==0) 
        {
          printf("\nDNA sequences detected. DNA GTR model will be used and optimisation of parameters carried out.\n");
        }
      if(datatype==1) 
        {
          printf("\nProtein sequences detected. Best Protein model will be determined.\n");
        }
      }
    }
    fclose(nexusfile);
  /* End test for datatype */

  /* initiate PAUP pipe */

  if(setup_pipe () == TRUE)
    {
    printf("All temporary files for this run will be annotated with the PID of this process (%d)\n", pid);
    interval1 = time(NULL);



    printf("\nExecuting file %s", infile);
    sprintf(command, "execute %s; exclude uninf;", infile);
    send_command (command);
    checkcount=0; while((response = check_for_output("Data matrix has")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    if(response == 3)
      {
      fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
      exit(0);
      }

    token = strtok(output, " ");
    for(i=0; i<3; i++) token = strtok(NULL, " ");
    number_taxa = atoi(token);
    for(i=0; i<2; i++) token = strtok(NULL, " ");
    alignment_length = atoi(token);
    printf("\n\tNumber of taxa = %d Number of Characters = %d\n", number_taxa, alignment_length);

    /* get number of characters excluded that were not parsimony infomrative */
    response = check_for_output("characters excluded");
    token = strtok(output, " ");  /* get past the inital spaces on this line */
    /* get the number of excluded characters */
    excludedchars=atoi(token);
    printf("\tNumber of parsimony uninformative characters = %d\n", excludedchars);
    printf("\tNumber of Characters remaining in alignment = %d\n\n", alignment_length-excludedchars);


    /* Check to see if the nexus file constained a tree */
    if(buildflag == FALSE)
      { 
      if((response = check_for_output("read from TREES block")) == 1) /* indicates that a trees block was provided with the nexus file */
        {
        num_trees = atoi(strtok(output, " "));
        printf("\t%d tree(s) provided in the nexus file, the first will be used for the likelihood decay calculations\n\n", num_trees);
        }
      }
    if(datatype==0) send_command ("set criterion=like; lset NST=6 lCollapse=no;");
    if(datatype==1) send_command ("set criterion=like; lset aaFreq=empirical aaRMatrix=JTT lCollapse=no;");

    /* generate arrays to hold the sites likelihood later */
    site_likelihoods = malloc((number_taxa-2)*sizeof(double *));
    for(i=0; i<(number_taxa-2); i++) site_likelihoods[i]=malloc((alignment_length-excludedchars)*sizeof(double));


    if(num_trees == 0)
      {  
      likelihood = build_starting_tree (likelihood);

      interval2 = time(NULL);
      print_time(difftime(interval2, interval1));
      }
      interval2 = time(NULL);

    likelihood = re_estimate_parameters (likelihood);

    interval3 = time(NULL);
    print_time(difftime(interval3, interval2));

    /* now read tree from PAUP file, capturing the name definitions */
    treefilename = malloc(1000*sizeof(char));
    treefilename[0] = '\0';
 /*   checkcount =0; while(!check_for_output("saved to file")) { if(checkcount == 0) { printf("\n\tWaiting for PAUP .."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);} */
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
    if(norevcon == FALSE) /* as long as we haven;t indeicated that we wish to skip the reverse constraint analysis, using the "-n" option */
      {
      /* Next identify each split and for each, define a constraint, and fine the likelihood of the best tree that does not contain that split */
        printf("\nExtracting constraints defined by best tree and carry out reverse constraint analysis for each\n");

        interval2 = time(NULL);
        test_reverse_constraints(newtree);

        if(listcons==FALSE)
          {   
          interval3 = time(NULL);
          print_time(difftime(interval3, interval2));
     

          sprintf(command, "%s.likelihood.decay.tre", infile);
          outputtree = fopen(command, "w");
          fprintf(outputtree, "%s [%g]\n", newtree, likelihood);
          fclose (outputtree);

          printf("\n\n\nresulting weighted tree:\n%s [%g]\n\nTree with labels saved to file %s\n\n", newtree, likelihood, command);
          }
     /* Clean up and exit */
      free(treefilename);

      /*unroottree(string); */

      printf("\tFinished revered constraints analysis.\n");
      print_time(difftime(time(NULL), interval2));

      
      sprintf(str, "%s.sitelike.txt", infile);
      file=fopen(str, "w");
      fprintf(file, "Best_tree\t");
      for(i=0; i<(number_taxa-3); i++) fprintf(file, "constraint_tree %d\t", i);
      fprintf(file, "\n");

      for(j=0; j<alignment_length-excludedchars; j++)
        { 
        for(i=0; i<(number_taxa-2); i++)
          {
          fprintf(file, "%f\t", site_likelihoods[i][j]);
          }
        fprintf(file, "\n");
        }
      fclose(file);


      }

    if(bootreps > 0) /* carry out a boostrap analysis - NOTE: this is only for illustrative purposes and the resulting topology is NOT used for the reverse constraint analysis (the "best" likelihood tree os used for the REverse constraint analysis) */
      { 
      interval2 = time(NULL);
      do_bootstrap_analysis (bootreps);
      interval3 = time(NULL);
      print_time(difftime(interval3, interval2));
      }

    error = close_pipe();
    printf("Finished analysis, Overall ");
    print_time(difftime(interval3, interval1));
    }
  else
    {
    error = TRUE;
    }




  if(deltmp==TRUE)
    { 
    sprintf(str, "rm -f Paup_output_%s.txt", pid_string);
    system(str);
    sprintf(str, "rm -f paup_%s.altnexus.tre", pid_string);
    system(str);   
    sprintf(str, "rm -f paup_%s.tre", pid_string);
    system(str);
    sprintf(str, "rm -f paup_%s_log.txt", pid_string);
    system(str);
    sprintf(str, "rm -f paupblock_%s.command", pid_string);
    system(str);
    sprintf(str, "rm -f site_like_scores_%s.txt", pid_string);
    system(str);
    sprintf(str, "rm -f site_like_scores_%s.txt", pid_string);
    system(str);
    }
  

 
  if(gtflag == TRUE) fclose(guidetreefile);

  for(i=0; i<(number_taxa-2); i++) free(site_likelihoods[i]);
  
  free(site_likelihoods);
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
  free(newtree);
  free(guidetree);
  free(str);

return(error);
}







int check_for_output ( char *search_string )
  {
  int found=FALSE, error=FALSE, i=0;


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
    while(!feof(logfile) && !found && !error) {
          output[i]= getc(logfile);
          i++;
          if(output[i-1] == '\n' || output[i-1] == '\r')
          {
          output[i] = '\0';
          if(strstr(output, search_string ) != '\0')
              found = 1;
          if(strstr(output, "Error" ) != '\0')
              error=TRUE;
          i=0;
          }
      }
    }
  fclose(logfile);

  if(error == TRUE) found =3;
  /*** End looking for output */
  return(found);
  }



void send_command(char * sent_command)
  {
  char com[1000], paupblockfile[10000];
  FILE *paupblock = NULL;

  com[0] = '\0'; paupblockfile[0] = '\0';

  /* delete logfile, so that outputs from previous commands do not get mixed up with the outputs of the command about to be run */
  remove(logfilename);
  strcpy(com, "touch ");
  strcat(com, logfilename);
  system(com);

  sprintf(paupblockfile, "paupblock_%d.command", pid);
  paupblock=fopen(paupblockfile, "w");
  fprintf(paupblock, "#NEXUS\nbegin paup;\n%s\nend;\n", sent_command);
  fflush(paupblock);
  fclose(paupblock);
  if(print_command == TRUE) fprintf(stderr, "\n\t\t[Paup command: %s]\n",sent_command );
  fflush(stdout);

 /* printf ("Command sent: %s\n", sent_command); */
  fprintf (paup_pipe, "log file=%s start=yes replace=yes flushLog=yes; \n", logfilename);
  fflush(paup_pipe);


  fprintf (paup_pipe, "exe %s;\n", paupblockfile);
  fflush(paup_pipe);
  fprintf (paup_pipe, "log stop=yes;\n");
  fprintf (paup_pipe, ";");
  fflush(paup_pipe);
  sleep(2);
  }


double re_estimate_parameters (double likelihood)
  {
  double new_likelihood = 0, test_likelihood = likelihood, prev_likelihood = likelihood ;
  int i=0, iteration=1, best_aafeq=0, best_aaRMatrix=0, best_rates=0, stop=FALSE;
  char * token, *comm;
  size_t t = 0, q=0;
  char aafreq[][10] = { "empirical", "equal" }, aaRMatrix[][10] = { "JTT", "JTTPAML", "PAM", "MTrev24", "WAG", "LG" }, rates[][30] = { "equal", "gamma shape=estimate" };

  comm = malloc(10000*sizeof(char));
  comm[0] = '\0'; 

      if(num_trees > 0)
        { 
          printf("\n\tDeterming likelihood for given tree\n");
          send_command("lscores 1;");
          checkcount=0; while((response = check_for_output("-ln L")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
          if(response == 3)
            {
            fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
            exit(0);
            }
          token = strtok(output, " ");
          for(i=0; i<2; i++) token = strtok(NULL, " ");
          likelihood=atof(token);
          printf( "\n\tBest -ln L: %g\n", likelihood );
          test_likelihood = likelihood;
          prev_likelihood = likelihood;
        }


      /* First determine if we need to use a gamma distribution. */
      printf ("Begin testing for optimal site rates:\n");

      for( t = 0; t < sizeof(rates) / sizeof(rates[0]); t++)
        {
        sprintf(comm, "lset rates=%s; lscores 1;", rates[t]);  

        printf("\n\tCalculating likelihood with rates=%s", rates[t]); fflush(stdout);
        send_command (comm);
        checkcount=0; while((response = check_for_output("-ln L")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
        if(response == 3)
          {
          fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
          exit(0);
          }
        token = strtok(output, " ");
        for(i=0; i<2; i++) token = strtok(NULL, " ");
        test_likelihood=atof(token);
        printf( "\n\tBest -ln L: %g\n", test_likelihood );
        if((fabs(round(test_likelihood)) - fabs(round(prev_likelihood))) < 0)
          {
          prev_likelihood = test_likelihood;
          best_rates = t;
          }
        }
      if(best_rates==1)
        sprintf(comm, "lset rates=gamma shape=previous;");
      else
        sprintf(comm, "lset rates=%s;", rates[best_rates]);
      send_command (comm);
      printf("\n\tBest rates model found: %s\n\n", rates[best_rates]);


      likelihood=prev_likelihood;

      if(datatype==1) /* Protein sequences */
        {

          /* Estimate the best aafreq ( from "empirical|equal|estimate|") */
        printf ("Begin testing for optimal AA Rmatrix and frequencies:\n");
        while((fabs(round(new_likelihood)) - fabs(round(likelihood))) < 0 && stop == FALSE)
          {
          if(iteration>1) 
              {
                printf("\tTree with better likelihood found using optimised models, reestimation best models with this tree\n\n");
                likelihood = new_likelihood;
              }
          prev_likelihood=likelihood;
          for( t = 0; t < sizeof(aafreq) / sizeof(aafreq[0]); t++)
            {
            for( q = 0; q < sizeof(aaRMatrix) / sizeof(aaRMatrix[0]); q++)
              {
           /*   sprintf(comm, "lset rates=gamma shape=estimate aafreq=%s aaRMatrix=%s; lscores;", aafreq[t], aaRMatrix[q]);  */
              sprintf(comm, "lset aafreq=%s aaRMatrix=%s; lscores 1;", aafreq[t], aaRMatrix[q]);  

              printf("\n\tCalculating likelihood with aafreq=%s and aaRMatrix=%s", aafreq[t], aaRMatrix[q]); fflush(stdout);
              send_command (comm);
              checkcount=0; while((response = check_for_output("-ln L")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
              if(response == 3)
                {
                fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
                exit(0);
                }
              token = strtok(output, " ");
              for(i=0; i<2; i++) token = strtok(NULL, " ");
              test_likelihood=atof(token);
              printf( "\n\t-ln L: %g\n", test_likelihood );
              if(fabs(round(test_likelihood)) - fabs(round(prev_likelihood)) < 0)
                {
                prev_likelihood = test_likelihood;
                best_aafeq = t; best_aaRMatrix=q;
                }
              }
            }

          printf("\nBest AA model found: aafreq=%s, aaRMatrix=%s\n", aafreq[best_aafeq], aaRMatrix[best_aaRMatrix]);
          printf("\tRe-estimating parameters with best model");
          if(strcmp(aafreq[best_aafeq], "estimate")== 0)
            sprintf(comm, "lset aafreq=%s aaRMatrix=%s; lscores 1; lset aafreq=previous;", aafreq[best_aafeq], aaRMatrix[best_aaRMatrix]);
           /* sprintf(comm, "lset rates=gamma shape=estimate aafreq=%s aaRMatrix=%s; lscores; lset shape=previous aafreq=previous;", aafreq[best_aafeq], aaRMatrix[best_aaRMatrix]); */
          else
            sprintf(comm, "lset aafreq=%s aaRMatrix=%s; lscores 1;", aafreq[best_aafeq], aaRMatrix[best_aaRMatrix]);
           /* sprintf(comm, "lset rates=gamma shape=estimate aafreq=%s aaRMatrix=%s; lscores; lset shape=previous;", aafreq[best_aafeq], aaRMatrix[best_aaRMatrix]); */

          send_command(comm);
          checkcount=0; while((response = check_for_output("-ln L")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
          if(response == 3)
            {
            fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
            exit(0);
            }
          token = strtok(output, " ");
          for(i=0; i<2; i++) token = strtok(NULL, " ");
          new_likelihood=atof(token);
          printf( "\n\tBest -ln L: %g\n", new_likelihood );

          if((fabs((round(new_likelihood))) - fabs(round(likelihood))) < 0 && num_trees ==0)
            { 
            printf("\nUsing predicted aafreq and searching for best tree");
            sprintf(comm, "hs;");
            send_command(comm);
            checkcount=0; while((response = check_for_output("Score of best tree(s) found")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
            if(response == 3)
              {
              fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
              exit(0);
              }
            token = strtok(output, " ");
            for(i=0; i<6; i++) token = strtok(NULL, " ");
            new_likelihood=atof(token);
            printf( "\n\tBest -ln L: %g\n", new_likelihood );
            }
          if(num_trees > 0) stop=TRUE;
          iteration++;
          }
        if(num_trees == 0)
          printf("\tDifference in -ln L less than 1 unit, stopping search\n");
        else
          printf("\tCalculation of likelihood of given tree complete\n");
        }


      if(datatype==0) /* DNA sequences */
        {

         while((fabs((round(new_likelihood))) - fabs(round(likelihood))) < 0 && stop == FALSE)
            {
            if(iteration>1) 
              {
                printf("\tTree with changed likelihood found, starting reestimation of paramenters with this tree\n\n");
                likelihood = new_likelihood;
              }

            printf("Iteration %d: Estimating parameters given tree and recalculating likelihood", iteration);

     
            send_command ("lset lCollapse=no Rmatrix=estimate Basefreq=Estimate; lscores 1; lset Rmatrix=previous Basefreq=previous;");
           /* send_command ("lset Rmatrix=estimate Basefreq=Estimate rates=gamma shape=estimate; lscores; lset Rmatrix=previous Basefreq=previous shape=previous;"); */

            checkcount=0; while((response = check_for_output("-ln L")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
            if(response == 3)
              {
              fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
              exit(0);
              }


            token = strtok(output, " ");
            for(i=0; i<2; i++) token = strtok(NULL, " ");
            new_likelihood=atof(token);
            printf( "\n\tBest -ln L: %g\n", new_likelihood );
        
                     
            if((fabs((round(new_likelihood))) - fabs(round(likelihood))) < 0 && num_trees == 0)
              {
              likelihood = new_likelihood;

              printf("Iteration %d: Changed likelihood detected. Beginning tree search", iteration);
              send_command ("hs;");
              checkcount=0; while((response = check_for_output("Score of best tree(s) found")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
              if(response == 3)
                {
                fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
                exit(0);
                }

              token = strtok(output, " ");
              for(i=0; i<6; i++) token = strtok(NULL, " ");
              new_likelihood=atof(token);
              printf( "\n\tBest -ln L: %g\n", new_likelihood );
              }
            if(num_trees > 0) stop=TRUE;
            iteration++;
            }
        if(num_trees == 0)
          printf("\tDifference in -ln L less than 1 unit, stopping search\n");
        else
          printf("\tCalculation of likelihood of given tree complete\n");
        }

    /* Calculate all the site likelihoods for the best tree found and read into memory */
    printf("\tCalculating site likelihoods of best tree");
    sprintf(comm, "lscores 1 /sitelikes=yes scoreFile=site_like_scores_%d.txt replace=yes display=no;", pid);
    send_command (comm);
    checkcount=0; while((response = check_for_output("Tree scores (and parameter estimates) written to file:")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    if(response == 3)
      {
      fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
      exit(0);
      }
    read_sitelike_file(-1);

    sprintf(comm, "savetrees file=paup_%d.tre brLens=yes;", pid);
    send_command(comm);
    sprintf(comm, "savetrees file=%s.bestMLtree.tre brLens=yes format=altnexus;", infile);
    send_command(comm); 
    printf("\nBest tree found after %d iterations with -ln L %g saved to file %s.bestMLtree.tre;\n\n", iteration-1,  new_likelihood, infile); 

    free(comm);
    return(new_likelihood);
  }
 

  double build_starting_tree (double likelihood)
    {
    char * token;
    int i;

    printf("Generating starting ML tree");
    
    send_command("hs;");
    checkcount=0; while((response = check_for_output("Score of best tree(s) found")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    if(response == 3)
      {
      fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
      exit(0);
      }

    token = strtok(output, " ");
    for(i=0; i<6; i++) token = strtok(NULL, " ");
    likelihood=atof(token);
    printf( "\n\tBest -ln L: %g\n\n", likelihood );

    return(likelihood);
    }

 void do_bootstrap_analysis (int reps)
    {
    char *comm=NULL;
    int i;

    comm = malloc(10000*sizeof(char));
    comm[0] = '\0'; 

    printf("Beginning ML Boostrap analysis with %d replicates\n", reps);

    sprintf(comm, "boot nrep=%d treefile=%s.boottrees.tre format=altnexus / enforce=no; exe %s.boottrees.tre; contree /LE50=yes majRule=yes strict=no treefile=%s.contree.tre;",reps, infile, infile, infile);
    send_command(comm );  

    checkcount=0; while((response = check_for_output("Consensus tree(s) written to treefile:")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    if(response == 3) 
      {
      fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
      exit(0);
      }
    printf("\n\tResulting boostrapped phylogeny written to file: %s.boottrees.tre\n", infile);
    printf("\tMajority-Rule consensus of boostrapped trees (with compatible minor components) written to file: %s.contree.tre\n", infile);

    free(comm);
    } 


 int setup_pipe (void)
 { 
  int success = TRUE;
  char start_comm[10000];

  start_comm[0] = '\0';
  sprintf(start_comm, "paup -n > Paup_output_%d.txt 2> /dev/null", pid);
  /* Open the pipe */
  paup_pipe = popen (start_comm, "w");

  /* Check that pipes are non-null, therefore open */
  if ((!paup_pipe))
    {
      fprintf (stderr,
               "Paup pipe failed.\n");
      success = FALSE;
    }
  else
    {  

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

    /* Stop Paup from printing output unnecessarily to the screen */
    send_command ("set increase=auto warnReset=no warnTree=no warnTSave=no warnBlkName=no warnRoot=no  warnRedef=no errorStop=no errorBeep=no queryBeep=no flushLog=yes taxLabels=full;");
 
    }

  return(success);
 }


int close_pipe ()
  {
  int success = TRUE;
   /* Quit Paup */
    send_command("quit;");

    /* Close paup_pipe, cehcking for errors */
    if (pclose (paup_pipe) != 0)
      {
        fprintf (stderr,
                 "Could not close 'paup', or other error.\n");
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
                          fprintf(stderr, "Error: taxa %s is not in the translation table\n", single);
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
    int i=0, k=0, l=0, j=0, x=0, y=0, z=0, q=0, splitnumber=0, constraint_num=0, numsites_supporting_constraint=0, numsites_supporting_neither=0;
    char taxaname[1000], constraint[1000000], weight[10000], tree[1000000], tmptree[1000000];
    double sumlikesuppbest=0, sumlikesuppconstr = 0;


    taxaname[0] = '\0';
    constraint[0] = '\0'; j=0;
    weight[0] = '\0';
    tmptree[0] = '\0';
    tree[0]='\0';

    total_constraints=0;

    strcpy(tree, translated_tree);
     while(tree[i] != ';')
      {
        if(tree[i] == '(') total_constraints++;
        i++;
      }
      total_constraints--;
      total_constraints--;
    if(constart < 1) constart = 0;
    if(constart > total_constraints) constart = total_constraints;
    if(conend == -1) conend = total_constraints;
    if(conend > total_constraints) conend = total_constraints;

    resultingtrees = malloc((total_constraints+2)*sizeof(char *));
    for(i=0; i<total_constraints+2; i++)
      {
      resultingtrees[i]=malloc(1000000*sizeof(char));
      resultingtrees[i][0]= '\0';
      }
    /* copy best tree to resultingtrees[0] */
    strcpy(resultingtrees[0], tree);


    i=0;
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
          
          if(constraint_num >= constart && constraint_num <= conend)
            /* Calculate the likeliood of the best tree that does not contain this constraint */
          {
            if(listcons==FALSE)
              {
              new_likelihood = build_starting_constraint_tree (likelihood, constraint, constraint_num);

             /* new_likelihood = re_estimate_constraint_parameters (new_likelihood, constraint_num); */ /* It is likely unnecessary to re-estimate all the parameters for the constraint trees when everything has been optimised already for the "best" tree */
              numsites_supporting_constraint =0; numsites_supporting_neither=0;
              sumlikesuppbest =0; sumlikesuppconstr=0;
             for(x=0; x<alignment_length-excludedchars; x++) 
              {  
              if((site_likelihoods[0][x] - site_likelihoods[constraint_num+1][x]) > 0 )
                {
                  numsites_supporting_constraint++;
                  sumlikesuppconstr=sumlikesuppconstr+(site_likelihoods[0][x] - site_likelihoods[constraint_num+1][x]);
                }
              if((site_likelihoods[0][x] - site_likelihoods[constraint_num+1][x]) < 0 )
                {
                  sumlikesuppbest=sumlikesuppbest+(site_likelihoods[0][x] - site_likelihoods[constraint_num+1][x]);
                }
     
              if((site_likelihoods[0][x] - site_likelihoods[constraint_num+1][x]) == 0 ) numsites_supporting_neither++;
              }    

              printf("\n\tBest reverse Constraint -ln L = %g\n\tDifference = %g (composed of %g supporintg unconstrained tree and %g supporting constrained tree)\tProportion likelihood decay  = %g\n", new_likelihood, new_likelihood-likelihood, fabs(sumlikesuppbest), fabs(sumlikesuppconstr), (fabs(sumlikesuppbest)/((fabs(sumlikesuppconstr)+(fabs(sumlikesuppbest))))));
             /* sprintf(weight, "%g/%g", new_likelihood- likelihood, (new_likelihood- likelihood)/(meanRand_likelihood-likelihood)); */
              sprintf(weight, "%d/%g/%g/%g/%g", constraint_num, new_likelihood- likelihood, fabs(sumlikesuppbest), fabs(sumlikesuppconstr) , (fabs(sumlikesuppbest)/((fabs(sumlikesuppconstr)+(fabs(sumlikesuppbest))))));


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
              }  
            else
              {
              printf("Constraint %d = %s\n", constraint_num, constraint);

              }
          }


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
    printf("\t\n====================\nSetting constraint number %d/%d: %s\n====================", constraint_num, total_constraints, constraint);
    sprintf(comm, "constraints constraint_%d = (%s);", constraint_num, constraint);
    send_command(comm );

    printf("\n\tGenerating ML reverse constraint tree");
    sprintf(comm, "hs constraints=constraint_%d enforce=yes converse=yes;", constraint_num);
    send_command (comm);
 
    checkcount=0; while((response = check_for_output("Score of best tree(s) found")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    if(response == 3)
      {
      fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);
      exit(0);
      }
    token = strtok(output, " ");
    for(i=0; i<6; i++) token = strtok(NULL, " ");
    likelihood=atof(token);
    printf( "\n\tBest -ln L: %g\n", likelihood );


    /* Calculate all the site likelihoods for the best tree found and read into memory */
    printf("\tCalculating site likelihoods of reverse constraints tree\n");
    sprintf(comm, "lscores 1 /sitelikes=yes scoreFile=site_like_scores_%d.txt replace=yes display=no;", pid);
    send_command (comm);
    checkcount=0; while((response = check_for_output("Tree scores (and parameter estimates) written to file:")) != 1 && response != 3 )  { if(checkcount == 0) { printf("\n\tWaiting for PAUP ..."); fflush(stdout); checkcount++;} sleep(5); printf("."); fflush(stdout);}
    if(response == 3)
      {
      fprintf(stderr, "Error signal detected from Paup. Now quitting. Please check the file \"Paup_output_%d.txt\" to determine the error message\n", pid);      
      exit(0);
      }
    read_sitelike_file(constraint_num);

    if(constraint_num == 0)
        {
        sprintf(sys, "rm -f %s.constraint.tre;",infile );
        system(sys);
        }
    /* append the best reverse constraint tree to the output file */
    sprintf(sys, "echo \"[constraint tree %d. -ln L %g]\" >> %s.constraint.tre", constraint_num, likelihood, infile);
    system(sys);
    sprintf(comm, "savetrees file=%s.constraint.tre brLens=yes format=altNexus append=yes;",infile );
    send_command (comm);
    
    

    free(comm);
    return(likelihood);
    }


  void print_time(double diff)
    {
        if(diff > 0)
            {
            printf("Time taken:"); 
            if(diff > 60)
                {
                if(diff > 3600)
                    {
                    if(diff > 86400)
                        {
                        printf(" %0.0lf Day", diff/86400);
                        if(diff/86400 > 1) printf("s");
                        diff = diff- (86400 * ((int)diff/86400));
                        }
                    printf(" %0.0lf Hour", diff/3600);
                    if(diff/3600 > 1) printf("s");
                    diff = diff - (3600 * ((int)diff/3600));
                    }
                printf(" %0.0lf Minute", diff/60);
                if(diff/60 > 1) printf("s");
                diff = diff - (60 * ((int)diff/60));
                }
            printf(" %0.0lf second", diff);
            if(diff > 1) printf("s");
            printf("\n\n");
            }
    }


void read_sitelike_file(int component) /* this will be passed "-1" to indicate that these are the site likelihhods for the "best" tree */
  {

  FILE *sitelike_file = NULL;
  char filename[10000], str1[100], str2[100], str3[100], str4[100], c;
  int line=0, part=0, numlineends=0;

  if(component == -1) part =0;
  else part=component+1;

  sprintf(filename, "site_like_scores_%d.txt", pid);
  sitelike_file = fopen(filename, "r");
   while(!feof(sitelike_file) && numlineends <2)
    {    
      if((c=getc(sitelike_file)) == '\n' || c == '\r') numlineends++;
    }

  while(!feof(sitelike_file))
    {
    fscanf(sitelike_file, "%s\t%s\n", str1, str2 );
      site_likelihoods[part][line] = atof(str2);
    line++;
    }
  fclose(sitelike_file);
  }
