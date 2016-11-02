/* Author: Ram Samudrala (me@ram.org)
 * March 1, 1996. 
 * modified Hong Hung Sept 2006 2012
 */

#include "lite.h"

/******************************************************************/

int check_atomp(atom *atomp, const char routine_name[])
{
  if (atomp == ((atom *) NULL))
    {
      fprintf(stderr, "%s(): atomp is NULL when it shouldn't be!\n", routine_name);
      exit(FALSE);
    }
  return TRUE;
}

/******************************************************************/

int check_eof(int status, const char routine_name[])
{
  if (status == EOF)
    fprintf(stderr, "check_eof(): EOF returned in %s()!\n", routine_name);
  return TRUE;
}

/******************************************************************/

int check_malloc(void *pointer, const char pointer_string[], const char routine_name[])
{
  if (pointer == ((void *) NULL))
    {
      fprintf(stderr, "%s(): unable to allocate memory for %s!\n", routine_name, pointer_string);
      exit(FALSE);
    }
  return TRUE;
}

/******************************************************************/

int check_maximum_value(int value, int maximum_value, const char routine_name[])
{
  if (value >= maximum_value)
    {
      fprintf(stderr, "%s(): value, %d, exceeded maximum_value (%d)!\n", 
	      routine_name, value, maximum_value);
      exit(FALSE);
    }
  return TRUE;
}

/******************************************************************/

int check_null(void *pointer, const char routine_name[])
{
  if (pointer == ((void *) NULL))
    {
      fprintf(stderr, "%s(): NULL encountered unexpectedly!\n", routine_name);
      exit(FALSE);
    }
  return TRUE;
}

/**
 * @brief Open binary trajectory in AMBER7 binpos format for reading
 *
 * @param file_pp : pointer to pointer to the file
 * @param filename : name of the trajectory file
 * @param size_p: pointer to variable sizing up all the coordinates
 * @param nat_p : pointer to variable number of atoms in a frame
 * @param nmodels_p: pointer to variable number of models in list of names
 */
int open_binpos_read(FILE **file_pp, const char filename[], size_t *size_p, int *nat_p, int *nmodels_p)
{
  FILE *file_p;
  char magicchar[5];
  size_t point;
  int nat;
  char lenbuf[4];
  char tmpc;
  int er=0;
  size_t size=0;
  int nmodels=*nmodels_p;
  size_t record_size=0;
  size_t magic_size = (size_t)(4*sizeof(char));

  /* Open file for reading */
  file_p = fopen(filename, "rb");
  if (!file_p)
  {
    fprintf(stderr, "Could not open file '%s' for reading.\n", filename);
    exit(FALSE);
  }

  /* Compute the memory occupied by all binpos records */
  fseek(file_p , 0 , SEEK_END);
  size = (size_t)(ftell(file_p) - magic_size);
  *size_p = size;
  rewind (file_p);

  /*  Read magic number */
  fread(magicchar, sizeof(char), 4, file_p);
  magicchar[4]= '\0' ;
  if(strcmp(magicchar, "fxyz")!=0)
  {
    fprintf(stderr, "not a binpos amber coordinate file\n");
    exit(FALSE);
  }

  /* Read number of atoms */
  fread(&nat, sizeof(int), 1, file_p);
  fprintf(stderr, "Number of atoms is %d\n", nat);
  /* Check for endianness here */
  if(nat>1000000000){
    fprintf(stderr, "File '%s' appears to be other-endian.\n", filename);
    exit(FALSE);
  }
  *nat_p = nat;
  *file_pp = file_p;

  /* Check number of models in the list file is same as number of frames in the trajectory file */
  record_size = sizeof(int)+3*nat*sizeof(float);
  if(nmodels != (int)(size/record_size)){
    fprintf(stderr,"List of names contains %d names but trajectory contains %d frames\n",
            nmodels, size/record_size);
    exit(FALSE);
  }
  /* check size is an exact multiple of the number of models in the list file */
  if(size != nmodels*record_size){
    fprintf(stderr,"non_integer number of models : size=%d but nmodels*record_size=%d\n",
            size, nmodels*record_size);
    exit(FALSE);
  }

  fprintf(stderr,"nat=%d nmodels=%d size=%d sizeof(float)=%d sizeof(int)=%d\n",
          nat, nmodels, size, sizeof(float), sizeof(int));
  return TRUE;
}

/**
 * @brief Read coordinates from AMBER7 binpos file to array
 *
 * @param file_p : pointer to the file
 * @param array_p : pointer to array storing coordinates
 * @param pdb_size: numer of coordinates per frame
 * @param nmodels : number of frames
 */
 void binpos_read(FILE *file_p, float *array_p, int pdb_size, int nmodels){
  size_t array_index = 0;
  size_t lapse = (size_t)(pdb_size);
  float *array_current_pointer = array_p;
  int nat;

  for(int i=1; i<=nmodels; i++) {
    /* Read coordinates for current frame */
    if (fread(array_current_pointer, sizeof(float), pdb_size, file_p) != lapse) {
      fprintf(stderr, "Failure reading coordinates from amber7 binary file.\n");
      exit(FALSE);
    }
    array_current_pointer += lapse;
    /* Read number of atoms */
    if (((fread(&nat, sizeof(int), 1, file_p)) != 1) && i<nmodels) {
      fprintf(stderr, "Failure reading number of atoms from amber7 binary file.\n");
      exit(FALSE);
    }
  }
}

/******************************************************************/

int open_file(FILE **fp, const char filename[], const char status[], const char routine_name[])
{
  char buf[20];

  if (strcmp(filename, STDIN_FILENAME) == 0)
    {
      *fp = stdin;
      return TRUE;
    }
  if (strcmp(filename, STDOUT_FILENAME) == 0)
    {
      *fp = stdout;
      return TRUE;
    }
  if (strcmp(filename, STDERR_FILENAME) == 0)
    {
      *fp = stderr;
      return TRUE;
    }

  switch (status[0])
    {
    case 'r':
      strcpy(buf, "1"
          "reading");
      break;
    case 'a':
      strcpy(buf, "appending");
      break;
    case 'w':
      strcpy(buf, "writing");
      break;
    default:
      fprintf(stderr, "open_file(): unknown status (%s) encountered!\n", status);
      break;
    }

  if ((*fp = fopen(filename, status)) == NULL)
    {
      fprintf(stderr, "%s(): couldn't open %s for %s!\n", routine_name, filename, buf);
      exit(FALSE);
    }
   else if(routine_name)
    fprintf(stderr, "%s(): opening %s for %s...\n", routine_name, filename, buf);

  return TRUE;
}

/******************************************************************/

int close_file(FILE **fp, const char filename[], const char routine_name[])
{
  if ((strcmp(filename, STDOUT_FILENAME) != 0) && (strcmp(filename, STDIN_FILENAME) != 0) &&
      (strcmp(filename, STDERR_FILENAME) != 0))
    {
     if(routine_name) fprintf(stderr, "%s(): closing %s.\n", routine_name, filename);
      fclose(*fp);
    }
  return TRUE;
}

/******************************************************************/

int myisnan_(int *nanflag, double *value)
{
  if (isnan(*value) == 0)
    *nanflag = 0;
  else
    *nanflag = 1;

  return TRUE;
}

/******************************************************************/
