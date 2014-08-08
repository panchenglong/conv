#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include "convcode.h"

namespace ratelesscode
{

// --------- private functions -------------

void ConvolutionalCode::init_code()
{
    int number_inputs = 1 << k;
    int number_ouputs = 1 << n;
    int number_memory = 0;

    for (int i = 0; i < memory.size(); i++) 
        number_memory += memory(i);

    assert(number_memory < 32);
    int number_states = 1 << number_memory; // NOTE! int maybe not enough for large memory

    next_state_table.set_size(number_states, number_inputs);
    output_table.set_size(number_states, number_inputs);

    itpp::ivec shift_register(k);

    for (int current_state = 0; current_state < number_states; current_state++) {
        for (int current_input = 0; current_input < number_inputs; current_input++) {
            output_table(current_state, current_input) = 0;
            state_to_regs(current_state, shift_register);

            append_input(current_input, shift_register);

            for (int j = 0; j < n; j++) {
                int temp = 0;
                for (int i = 0; i < k; i++) {
                    temp += hamming_weight(memory(i) + 1,
                                           shift_register(i) & gen_matrix(i, j)
                                           );
                }
                output_table(current_state, current_input) += itpp::pow2i(n - j - 1) * (temp % 2);
            }
            next_state_table(current_state, current_input) = next_state(shift_register);
        }
    }

    if (debug & PRINT_OUTPUT_TABLE) {
        std::cout << "output_table:\n" << output_table << std::endl;
    }
    if (debug & PRINT_NEXT_STATE_TABLE) {
        std::cout << "next_state_table" << next_state_table << std::endl;
    }
}

void ConvolutionalCode::state_to_regs(int state, itpp::ivec& regs)
{
    int i, mask = 0, shift = 0, tmp = state;

    for (i = k - 1; i >= 0; i--) {
        mask = (1 << memory(i)) - 1;
        tmp >>= shift;
        regs(i) = tmp & mask;
        shift += memory(i);
    }
}

void ConvolutionalCode::regs_to_state(const itpp::ivec& regs, int& state)
{
    int i, shift = 0, tmp;

    for (i = k - 1; i >= 0; i--) {
        tmp = regs(i) << shift;
        state |= tmp;
        shift += memory(i);
    }
}

void ConvolutionalCode::append_input(int input, itpp::ivec& regs)
{
    for (int i = 0; i < k; i++) 
        regs(i) |= (((input >> (k - i - 1)) & 1) << memory(i));
}

int ConvolutionalCode::next_state(itpp::ivec& regs)
{
    int state = 0;

    for (int i = 0; i < k; i++) 
        regs(i) >>= 1;
    regs_to_state(regs, state);

    return state;
}


void ConvolutionalCode::print_nodes_status(std::vector<node> &nodes)
{
    std::cout << "-------NODES INFO----------" << std::endl;
    for (int i = 0; i < nodes.size(); i++) {
        if (!nodes[i].explored) {
            break;
        }
        std::cout << "Node " << i << std::endl;
        std::cout << "*state = " << nodes[i].state << std::endl;
        std::cout << "*gamma = " << nodes[i].gamma << std::endl;
        std::cout << "*metrics = ";
        for (int j = 0; j < nodes[i].metrics.size(); j++) {
            std::cout << nodes[i].metrics[j] << " ";
        }
        std::cout << "\n";
        std::cout << "*tm = ";
        for (int j = 0; j < nodes[i].tm.size(); j++) {
            std::cout << "(" << nodes[i].tm[j].first << ", "
                << nodes[i].tm[j].second << ")";
        }
        std::cout << "\n*branch = " << nodes[i].branch << std::endl;
    }
}

// --------- public functions --------------

void ConvolutionalCode::encode(const itpp::bvec& uncoded_bits, itpp::bvec& coded_bits, bool tail)
{
    itpp::bvec new_uncoded_bits;

    if (tail) {
        // add zeros after the uncoded bits
        int max_m = itpp::max<int>(memory);
        itpp::bvec added_zeros(max_m * k);
        added_zeros.zeros();
        new_uncoded_bits = itpp::concat<itpp::bin>(uncoded_bits, added_zeros);
    } else 
        new_uncoded_bits = uncoded_bits;

    int size = new_uncoded_bits.size();
    coded_bits.set_size((size / k) * n);

    int state = 0;
    int input;
    for (int i = 0, j = 0; i < size; i += k, j += n) {
        input = bvec2int(new_uncoded_bits, i, k);
        itpp::bvec output = int2bvec(output_table(state, input), n);
        for (int k = 0; k < n; k++) 
            coded_bits(j + k) = output(k);
        state = next_state_table(state, input);
    }
}


// Fano decoding algoritm is used to implement the decoder.
int ConvolutionalCode::decode(const itpp::vec& received_signal, itpp::bvec& decoded_bits,
                              int delta, unsigned long max_cycles, bool tail)
{
    long t;             /* Threshold */
    long ngamma;

    int size = received_signal.size() / n;

    std::vector<node> nodes(size + 1);
    typedef std::vector<node>::iterator node_iter;

    node_iter start_node = nodes.begin();
    node_iter current_node;
    node_iter last_node = nodes.end() - 1;
    int max_m = itpp::max<int>(memory);
    node_iter tail_node = last_node - max_m;

    int number_outputs = 1 << n;
    int number_inputs = 1 << k;

    /* Compute all possible branch metrics for each symbol pair
     * This is the only place we actually look at the raw input symbols
     */
    int idx = 0;
    for (node_iter iter = start_node; iter < last_node; iter++) {
        if (debug & PRINT_METRICS) {
            std::cout << "metrics in  node " << iter - start_node << std::endl;
        }
        iter->explored = false;

        for (int i = 0; i < number_outputs; i++) {
            long m = 0, k = i;
            for (int j = n - 1; j >= 0; j--) {
                m += get_metric(k & 1, received_signal(idx + j));
                k >>= 1;
            }
            iter->metrics.push_back(m);

            if (debug & PRINT_METRICS) {
                std::cout << m << " ";
            }
        }
        idx += n;

        if (debug & PRINT_METRICS) {
            std::cout << "\n";
        }
        
    }

    current_node = start_node;
    current_node->state = 0;
    current_node->explored = true;

    /* Compute and sort branch metrics from root node */
    for (unsigned int i = 0; i < number_inputs; i++) {
        std::pair<long, unsigned int> p(current_node->metrics[output_table(current_node->state, i)], i);
        current_node->tm.push_back(p);
        
    }

    std::sort(current_node->tm.begin(), current_node->tm.end(),
              [](std::pair<long, unsigned int> p1, std::pair<long, unsigned int> p2){ return p1.first > p2.first; });

    current_node[1].state = next_state_table(current_node->state, current_node->tm[0].second);
    current_node->branch = 0; //Start with best branch

    unsigned long maxcycles = max_cycles * size;
    current_node->gamma = t = 0;

    /* Start the Fano decoder */
    unsigned long cnt;
    for (cnt = 1; cnt <= maxcycles; cnt++) {
        if (debug & PRINT_NODES_CONTENT) {
            std::cout << "ROUND = " << cnt << "Threshold = " << t << std::endl;
            print_nodes_status(nodes);
        }

        /* Look forward */
        ngamma = current_node->gamma + current_node->tm[current_node->branch].first;
        if (ngamma >= t) {
            /* Node is acceptable */
            if (current_node->gamma < t + delta) {
                /* First time we've visited this node;
                 * Tighten threshold.
                 */
                while (ngamma >= t + delta) t += delta;
            }
            /* Move forward */
            current_node[1].gamma = ngamma;

            if (++current_node == last_node) break;  /* Done! */
            current_node->explored = true;

            /* Compute and sort metrics, starting with the 
             * zero branch
             */

            if (current_node >= tail_node) {
                /* The tail must be all zeroes, so don't even
                 * bother computing the 1-branches there.
                 */
                std::pair<long, unsigned int> p(current_node->metrics[output_table(current_node->state, 0)], 0);
                current_node->tm.clear();
                current_node->tm.push_back(p);
                
            } else {
                current_node->tm.clear();
                for (int i = 0; i < number_inputs; i++) {
                    std::pair<long, unsigned int> p(current_node->metrics[output_table(current_node->state, i)], i);
                    current_node->tm.push_back(p);
                }

                std::sort(current_node->tm.begin(), current_node->tm.end(),
                          [](std::pair<long,unsigned int> p1, std::pair<long, unsigned int> p2) { return p1.first > p2.first; });
                current_node[1].state = next_state_table(current_node->state, current_node->tm[0].second);

            }
            current_node->branch = 0;   /* Start with best branch */
            continue;
        }
        /* Threshold violated, can't go forward */
        for (;;) {
            /* Look backward */
            if (current_node == start_node || current_node[-1].gamma < t) {
                /* Can't back up either.
                 * Relax threshold and and look
                 * forward again to better branch.
                 */
                t -= delta;
                if (current_node->branch != 0) {
                    current_node->branch = 0;
                }
                break;
            }
            /* Back up */
            if (--current_node < tail_node && current_node->branch != number_inputs - 1) {
                /* Search next best branch */
                current_node->branch++;
                break;
            } /* else keep looking back */
        }
    }

    /* Copy decoded data to user's buffer */
    if (cnt >= maxcycles) 
        return -1;  /* Decoder timed out */
    
    if (tail) {
        size -= max_m;
    }
    decoded_bits.set_size(size * k);

    for (int i = 0; i < size; i++) {
        int branch = nodes[i].branch;
        int input = nodes[i].tm[branch].second;
        for (int j = k - 1; j >= 0; j--) {
            decoded_bits(i * k + j) = itpp::bin(input & 1);
            input >>= 1;
        }
    }

    return 0;       /* Successful completion */
}

void ConvolutionalCode::set_decoder_env(int pamp, double pnoise, double pbias, int pscale, int pdebug)
{
    amp = pamp;
    noise = pnoise;
    bias = pbias;
    scale = pscale;
    debug = pdebug;
}

long ConvolutionalCode::get_metric(int snd, double received)
{
    double p0 = normal_distribution(amp * 1.0, noise, received);
    double p1 = normal_distribution(amp * -1.0, noise, received);

    double numerator;
    if (snd == 0) 
        numerator = p0;
    else
        numerator = p1;

    return std::floor((std::log2((2 * numerator)/(p0 + p1)) - bias) * scale  + 0.5);
}


// ---------- auxiliary functions-------------
int bvec2int(const itpp::bvec& b, int start, int len)
{
    int res = 0;
    for (int i = 0; i < len; i++) 
        res += itpp::pow2i(len - i - 1) * int(b[start + i]);
    return res;
}

itpp::bvec int2bvec(int number, int len)
{
    itpp::bvec res(len);
    int temp = number;
    for (int i = len - 1; i >= 0; i--) {
        res(i) = itpp::bin(temp & 1);
        temp >>= 1;
    }
    return res;
}


/*!
  Calculate the Hamming weight of the binary representation of in of size length
*/
int hamming_weight(int len, int in)
{
    int i, w = 0;
    for (i = 0; i < len; i++) {
        w += (in & (1 << i)) >> i;
    }
    return w;
}

}
