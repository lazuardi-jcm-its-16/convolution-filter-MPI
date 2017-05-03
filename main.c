#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "processing_img.h"
#include "send_wrappers.h"
#include "recv_wrappers.h"
#include "initializations.h"

MPI_Comm CARTESIAN_COMM;
extern MPI_Datatype mpi_block, mpi_block_img;

int* process_img_wrapper(int* block, int block_width, int block_height, int filter_rounds, int convergence_option,
                        int convergence_rounds, int rank, int img_width, int img_height);
unsigned char* convert_img(int* image, int img_height, int img_width);

/* Main arguments:
 * 1. the image to be processed path/<filename>
 * 2. the width of the image
 * 3. the height of the image
 * 4. number of times to apply the filter on the image
 * 5. is black and white (0/1)
 * 6. convergence (0/1)
 * 7. number of rounds to check for convergence
 */

int main(int argc, char *argv[])
{
    /*
    if ( argc < 7 )
    {
        printf("Wrong number of parameters given.\n");
        printf("Usage: %s <img_filename> <width> <height> <number of repetitions"
            "> <isBW image (0/1)> <withConvergence(0/1)> <number of rounds to "
            "check for convergence [optional]>\n", argv[0]);
        return 1;
    }
    int convergence_option = atoi(argv[6]);
    if ( convergence_option == 1 && argc != 8 )
    {
        printf("Wrong number of parameters given.\n");
        printf("Usage: %s <img_filename> <width> <height> <number of repetitions"
            "> <isBW image (0/1)> <withConvergence(0/1)> <number of rounds to "
            "check for convergence [optional]>\n", argv[0]);
        return 1;
    }
    char* filename = argv[1];
    int img_width = atoi(argv[2]);
    int img_height = atoi(argv[3]);
    int filter_rounds = atoi(argv[4]);
    int isBW = atoi(argv[5]);
    int convergence_rounds = -1;*/
    
    char* filename = "lena512.bmp";
    int img_width = 512;
    int img_height = 512;
    int filter_rounds = 1;
    int isBW = 1;
    int convergence_option = 0;
    int convergence_rounds = -1;
    
    if ( convergence_option == 1 )
        convergence_rounds = atoi(argv[7]);

    int numprocs, rank;

    // Read image and send it
    int *block, **rgb_block;
    if ( isBW )
        block = *(initalization_phase(filename, img_width, img_height, 1));
    else
        rgb_block = initalization_phase(filename, img_width, img_height, 0);

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    

    int block_width = img_width / sqrt(numprocs) + 2;
    int block_height = img_height / sqrt(numprocs) + 2;
    
    
    
    // Create new communicator (of Cartesian topology)
    int dims[2] = {(int) sqrt(numprocs), (int) sqrt(numprocs)};
    int cyclic[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, cyclic, 1, &CARTESIAN_COMM);

    int coords[2];
    memset(coords, '\0', sizeof(coords));
    if (CARTESIAN_COMM != MPI_COMM_NULL) {
        MPI_Cart_coords(CARTESIAN_COMM, rank, 2, coords);
        MPI_Cart_rank(CARTESIAN_COMM, coords, &rank);
    }
    else {
        printf("Could not properly set CARTESIAN_COMM_WORLD\n");
        MPI_Finalize();
        return EXIT_FAILURE;
    }
    MPI_Barrier(CARTESIAN_COMM);
    double start_time = MPI_Wtime();

    if ( isBW )
    {
        block = process_img(block, block_width, block_height, filter_rounds, convergence_option, convergence_rounds);
        int* image = NULL;
        if (rank != 0) {
            MPI_Gather(block + block_width + 1, 1, mpi_block, image, 1, mpi_block_img, 0, CARTESIAN_COMM);
        }
        else {
            image = malloc(img_width * img_height * sizeof(int));
            MPI_Gather(block + block_width + 1, 1, mpi_block, image, 1, mpi_block_img, 0, CARTESIAN_COMM);
            double end_time = MPI_Wtime();

            printf("%f\n", end_time - start_time);
            
            /*
            // write BW image as raw file
            FILE* output = fopen("output.raw", "wb");
            unsigned char* img_buffer = malloc(img_height * img_width *
                                            sizeof(unsigned char));

            int i;
            for (i = 0; i < img_height * img_width; i++) {
                img_buffer[i] = (unsigned char) image[i];
            }
            fwrite(img_buffer, sizeof(char), img_height * img_width, output);
            fclose(output);*/
            
            
            FILE *f;
            unsigned char *img_buffer = NULL;
            int filesize = 54 + 3*img_width*img_height;  //w is your image width, h is image height, both int
            if( img_buffer )
                free( img_buffer );
            img_buffer = (unsigned char *)malloc(3*img_width*img_height);
            memset(img_buffer,0,sizeof(img_buffer));
            
            int i;
            for (i = 0; i < img_height * img_width; i++) {
                img_buffer[i] = (unsigned char) image[i];
            }
            
            unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
            unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
            unsigned char bmppad[1] = {0};

            bmpfileheader[ 2] = (unsigned char)(filesize    );
            bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
            bmpfileheader[ 4] = (unsigned char)(filesize>>16);
            bmpfileheader[ 5] = (unsigned char)(filesize>>24);

            bmpinfoheader[ 4] = (unsigned char)(       img_width    );
            bmpinfoheader[ 5] = (unsigned char)(       img_width>> 8);
            bmpinfoheader[ 6] = (unsigned char)(       img_width>>16);
            bmpinfoheader[ 7] = (unsigned char)(       img_width>>24);
            bmpinfoheader[ 8] = (unsigned char)(       img_height    );
            bmpinfoheader[ 9] = (unsigned char)(       img_height>> 8);
            bmpinfoheader[10] = (unsigned char)(       img_height>>16);
            bmpinfoheader[11] = (unsigned char)(       img_height>>24);

            f = fopen("output.bmp","wb");
            fwrite(bmpfileheader,1,14,f);
            fwrite(bmpinfoheader,1,40,f);
            for(i=0; i<img_height; i++)
            {
                fwrite(img_buffer+(img_width*(img_height-i-1)),3,img_width,f);
                fwrite(bmppad,1,(4-(img_width)%4)%4,f);
            }
            fclose(f);
            
            
            free(img_buffer);
            free(image);
        }
        free(block);
    }
    else
    {
        int** rgb_img_buffer = malloc( 3 * sizeof(int*) );
        int i;
        for ( i = 0; i < 3; i++)
        {
            rgb_img_buffer[i] = process_img_wrapper(rgb_block[i], block_width, block_height, filter_rounds, convergence_option,
                                                    convergence_rounds, rank, img_width, img_height);
        }
        if (rank == 0)
        {
            int i, j;
            double end_time = MPI_Wtime();
            printf("%f\n", end_time - start_time);

            unsigned char** image_buffer = malloc( 3 * sizeof(unsigned char*));
            for ( i = 0; i < 3; i++)
            {
                image_buffer[i] = convert_img(rgb_img_buffer[i], img_height, img_width);
                free(rgb_img_buffer[i]);
            }
            /*
            FILE* output = fopen("output.raw", "wb");
            unsigned char* image = malloc( 3 * img_height * img_width * sizeof(unsigned char));
            for ( i = 0 , j = 0; i < img_height * img_width; i++ , j += 3 )
            {
                image[j] = image_buffer[0][i];
                image[j+1] = image_buffer[1][i];
                image[j+2] = image_buffer[2][i];
                
            }
            fwrite(image, sizeof(char), 3 * img_height * img_width, output);
            fclose(output);*/
            
            
            
            FILE *f;
            unsigned char *image = NULL;
            int filesize = 54 + 3*img_width*img_height;  //w is your image width, h is image height, both int
            if( image )
                free( image );
            image = (unsigned char *)malloc(3*img_width*img_height);
            memset(image,0,sizeof(image));

            for ( i = 0 , j = 0; i < img_height * img_width; i++ , j += 3 )
            {
                image[j] = image_buffer[0][i];
                image[j+1] = image_buffer[1][i];
                image[j+2] = image_buffer[2][i];
                
            }

            unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
            unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
            unsigned char bmppad[3] = {0,0,0};

            bmpfileheader[ 2] = (unsigned char)(filesize    );
            bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
            bmpfileheader[ 4] = (unsigned char)(filesize>>16);
            bmpfileheader[ 5] = (unsigned char)(filesize>>24);

            bmpinfoheader[ 4] = (unsigned char)(       img_width    );
            bmpinfoheader[ 5] = (unsigned char)(       img_width>> 8);
            bmpinfoheader[ 6] = (unsigned char)(       img_width>>16);
            bmpinfoheader[ 7] = (unsigned char)(       img_width>>24);
            bmpinfoheader[ 8] = (unsigned char)(       img_height    );
            bmpinfoheader[ 9] = (unsigned char)(       img_height>> 8);
            bmpinfoheader[10] = (unsigned char)(       img_height>>16);
            bmpinfoheader[11] = (unsigned char)(       img_height>>24);

            f = fopen("output.bmp","wb");
            fwrite(bmpfileheader,1,14,f);
            fwrite(bmpinfoheader,1,40,f);
            for(i=0; i<img_height; i++)
            {
                fwrite(image+(img_width*(img_height-i-1)*3),3,img_width,f);
                fwrite(bmppad,1,(4-(img_width*3)%4)%4,f);
            }
            fclose(f);
            
            for ( i = 0; i < 3; i++)
            {
                free(image_buffer[i]);
            }
            free(image_buffer);
            free(image);
            
        }
        free(rgb_img_buffer);
    }
    MPI_Finalize();
    return 0;
}

int* process_img_wrapper(int* block, int block_width, int block_height, int filter_rounds, int convergence_option,
                        int convergence_rounds, int rank, int img_width, int img_height)
{
    block = process_img(block, block_width, block_height, filter_rounds, convergence_option, convergence_rounds);
    int* image = NULL;
    if (rank != 0) {
        MPI_Gather(block + block_width + 1, 1, mpi_block, image, 1, mpi_block_img, 0, CARTESIAN_COMM);
        free(block);
        return NULL;
    }
    else {
        image = malloc(img_width * img_height * sizeof(int));
        MPI_Gather(block + block_width + 1, 1, mpi_block, image, 1, mpi_block_img, 0, CARTESIAN_COMM);
        free(block);
        return image;
    }
}

unsigned char* convert_img(int* image, int img_height, int img_width)
{
    unsigned char* img_buffer = malloc(img_height * img_width *
                                            sizeof(unsigned char));
    int i;
    for (i = 0; i < img_height * img_width; i++) {
        img_buffer[i] = (unsigned char) image[i];
    }
    return img_buffer;
}
