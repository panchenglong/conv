/**
 *  run this multiple time to see multiple different result. Or
 *  you can include in a loop here, all as your will.
 */

#include <iostream>
#include "convcode.h"

int main(int argc, char **argv)
{
    //std::cout << "------ test encoder ------" << std::endl;
//  itpp::ivec m =  "2 2";
//  itpp::imat g_matrix = "04 00 02; 00 04 03"; //Fig. 2. in the article <sequential decoding of convolutional codes>
    itpp::ivec m = "2";
    itpp::imat g_matrix = "07 05";
    ratelesscode::ConvolutionalCode code(m, g_matrix);

    double SNR = 10;
    double Eb = 1.0;
    double N0 = Eb / (std::pow(10, SNR / 10));
    itpp::vec symbols, rec;
    itpp::BPSK bpsk;
    itpp::BERC berc;  //The Bit Error Rate Counter class
                      //Count the number of errors

    itpp::RNG_randomize();

    itpp::bvec msg = itpp::randb(100);
    //itpp::bvec msg = "1 1 0 1";
    itpp::bvec coded_bits;
    code.encode(msg, coded_bits);

    std::cout << "coded bits:\n" << coded_bits << std::endl;


    //std::cout << "------- test fano decoder ------" << std::endl;


    //Do the BPSK modulation

    bpsk.modulate_bits(coded_bits, symbols);                                                        
                                                                                                    
    //std::cout << "modulated symbols = \n" << symbols << std::endl;                                  
                                                                                                    
    //Add the AWGN                                                                                  
    rec = std::sqrt(Eb) * symbols + std::sqrt(N0 / 2) * itpp::randn(coded_bits.size());             
    
    std::cout << "received signal = \n" << rec << std::endl;
                                                                                                    
    // PRINT_TABLES is the debug information , refer to the enum in the "convcode.h"                                                                                                     
    code.set_decoder_env(std::sqrt(Eb), std::sqrt(N0 / 2), code.get_rate(), 4, ratelesscode::PRINT_TABLES);                                                                                                                    
    itpp::bvec decoded_bits;                                                                        
    if (-1 == code.decode(rec, decoded_bits, 8, 100)) {
         std::cerr << "timed out" << std::endl;       
         return 1;
    }
                                                                                                    
    std::cout << "msg:\n" << msg << std::endl;                                                      
    std::cout << "decoded bits are (may be wrong)" << std::endl;                                    
    std::cout << decoded_bits << std::endl;                                                         
                                                                                                    
                                                                                                    
    berc.count(msg, decoded_bits);                                                                  
                                                                                                    
    //std::cout << "There were " << berc.get_errors() << " received bits in error." << std::endl;   
    //std::cout << "There were " << berc.get_corrects() << " correctly received bits." << std::endl;
    std::cout << "The error probability was " << berc.get_errorrate() << std::endl;                 

    //std::cout << "The theoretical error probability is " << 0.5 * erfc(1.0) << std::endl;
    //file << berc.get_errorrate() << std::endl;

    //file.close();

    return 0;
}
