#pragma once
#include <omp.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Common/block.h>
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtReceiver.h"
#include "libOTe/TwoChooseOne/Iknp/IknpOtExtSender.h"
#include "coproto/Socket/AsioSocket.h"
#include "libOTe/Base/BaseOT.h"
#include "cryptoTools/Common/BitVector.h"
#include "cryptoTools/Crypto/PRNG.h"
#include "coproto/Common/span.h"
#include "psi/Defines.h"
// #include "utils.h"

using namespace std;
using namespace osuCrypto;
using namespace coproto;
using namespace CmpFuzzyPSI;


template<typename T>
class NCN1OT {
public:
    // NCN1OT() = default;

    NCN1OT(Role role, int nums, int length, IknpOtExtSender *sender, IknpOtExtReceiver *recver, Socket &chl) {
        // this->party = party;
        // this->io = hsio[0];
        this->nums = nums;
        this->length = length;
        // cout << ot_bsize << endl;
        this->sender = sender;
        this->recver = recver;
        this->chl = chl;
        this->role = role;
        
        prg.SetSeed(sysRandomSeed());
    }

    ~NCN1OT() {
        // delete ferretcot;
    }

    void send(T** data) {
        int depth = log2(length);
        int numOTs = nums * depth;
        T **seeds = (T **)malloc(nums * sizeof(T *));
        T *str0 = new T[nums*depth], *str1 = new T[nums*depth];
        cout << "h" << nums*depth << endl;
        //#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < nums; ++i) {
            OblivSetup(length, seeds + i, str0+(i*depth), str1+(i*depth));
        }

        cout << "h" << endl;
        send_ot(str0, str1, nums*depth);

        cout << "h" << endl;
        //#pragma omp parallel for
        for (int i = 0; i < nums; ++i) {
            Expand(length, -1, *(seeds + i), data[i]);
        }
        cout << "h" << endl;

        delete[] str0;
        delete[] str1;
    }

    void recv(T **data, oc::span<u8> epsilon) {
        int depth = log2(length);
        int numOTs = nums * depth;
        T **seeds = (T **)malloc(nums * sizeof(T *));
        // bool *chosen_bit = new bool[nums*depth];
        BitVector chosen_bit(nums*depth);
        T *temp = new T[nums*depth];
        
        // PRNG prg;
        prg.get(epsilon);

        //#pragma omp parallel for
        for (int i = 0; i < nums; ++i) {
            OblivSetup(length, seeds + i, NULL, NULL);
            int tmp_i = i * depth;
            epsilon[i] = epsilon[i] % length;
            uint64_t x = epsilon[i];
            for (int j = depth-1; j >= 0; --j) {
                chosen_bit[tmp_i + j] = (x & 1) ^ 1;
                x >> 1;
            }
        }
        
        recv_ot(temp, chosen_bit, nums*depth);
        
        //#pragma omp parallel for
        for (int i = 0; i < nums; ++i) {
            memcpy(*(seeds+i), temp+i*depth, depth*sizeof(T));
            for (int j = 0; j < depth; ++j) {
            }
            Expand(length, epsilon[i], *(seeds + i), data[i]);
        }

        // delete[] chosen_bit;
        delete[] temp;
    }
private:
    // int party;
    Role role;
    Socket chl;
    IknpOtExtSender *sender;
    IknpOtExtReceiver *recver;
    PRNG prg;
    T *a, *b, *delta;
    // static int64_t ot_bsize;
    // MITCCRH<ot_bsize> mitccrh;

    int nums, length;
    void DLenPRG(T seed, T *out) {
        block temp = toBlock(0, seed);
        AES aes(temp);
        block res, pt = AllOneBlock;

        aes.ecbEncBlocks<1>(&AllOneBlock, &res);
        memcpy(out, &res, 2*sizeof(T));
    }

    void OblivSetup(uint64_t length, T **seeds, T *str0, T *str1) {
        int depth = int(log2(length));
        if (role == Role::Sender) {
            *seeds = (T *)malloc(1 * sizeof(T));

            // generate root
            // PRNG prg(sysRandomSeed());
            prg.get(*seeds, 1);
            // preparing data for the following ote
            // T str0[depth], str1[depth];
            memset(str0, 0x00, depth * sizeof(T));
            memset(str1, 0x00, depth * sizeof(T));

            PrepareCorrelation(0, depth, **seeds, str0, str1);

            // start ote
            // cout << "start ote" << endl;
            // send_ot(str0, str1, depth);
            // cout << "end ote" << endl;
        } else {
            *seeds = (T *)malloc(depth * sizeof(T));

            // parse input i to bool vector
            // bool chosen_bit[depth];
            // for (int i = depth - 1; i >= 0 ; i--) {
            //     chosen_bit[i] = (x & 1) ^ 1;
            //     x >>= 1;
            // }
            // cout << "start ote" << endl;
            // recv_ot(*seeds, chosen_bit, depth);
            // cout << "end ote" << endl;
        }
        // io->flush();
    }

    void SimpleExpand(int cur_depth, int depth, int index, T seed, T *v) {
        T next_seeds[2];
        DLenPRG(seed, next_seeds);

        if (cur_depth == depth - 1) {  
            v[index] = next_seeds[0];
            v[index + 1] = next_seeds[1];
            return;
        }

        int offset = 1 << (depth - cur_depth - 1);
        SimpleExpand(cur_depth + 1, depth, index, next_seeds[0], v);
        SimpleExpand(cur_depth + 1, depth, index + offset, next_seeds[1], v);
    }

    int CheckOnPath(T x, T y) {
        return x==y;
    }

    void PunctureExpand(int depth, uint64_t x, T *seeds, T *v) {
        int size = 0;
        int cur_depth = 1;

        bool bits[depth];

        T seed, agg_mask;
        T next_seeds[2];

        std::deque<T> q;

        for (int i = depth - 1; i >= 0 ; i--) {
            bits[i] = (x & 1) ^ 1;
            x >>= 1;
        }

        int path = bits[0];
        
        if (bits[0]) {
            q.push_back(0);
            q.push_back(*seeds);
        } else {
            q.push_back(*seeds);
            q.push_back(0);
        }

        while (cur_depth < depth) {
            size = q.size();
            agg_mask = 0;
            // process each layer
            for (int i = 0; i < size; i++) {
                seed = q.front();
                q.pop_front();
                if (!CheckOnPath(seed, 0)) { // not on-path node
                    // evaluate prg and push directly
                    DLenPRG(seed, next_seeds);
                    q.push_back(next_seeds[0]);
                    q.push_back(next_seeds[1]);
                    
                    // in addition, we need aggregate all seeds according to each layer's bit
                    agg_mask ^= next_seeds[bits[cur_depth]];
                } else {
                    // 0 plays a placeholder of punctured value
                    if (bits[cur_depth]) {
                        q.push_back(0);
                        q.push_back(seeds[cur_depth]);
                    } else {
                        q.push_back(seeds[cur_depth]);
                        q.push_back(0);
                    }
                }
            }
            // correct the seed 
            path = ((path ^ 1) << 1) ^ bits[cur_depth];
            q[path] ^= agg_mask;
            cur_depth++;
        }

        for (int i = 0; i < q.size(); i++) {
            v[i] = q[i];
        }
    }


    void Expand(uint64_t length, uint64_t x, T *seeds, T *v) {
        int depth = int(log2(length));
        // cout << "start expand" << v.size() << endl;
        if (role == Role::Sender) {
            SimpleExpand(0, depth, 0, *seeds, v);
        } else {
            PunctureExpand(depth, x, seeds, v);
        }
        // cout << "end expand" << endl;
    }

    void PrepareCorrelation(int cur_depth, int depth, T seed, T *str0, T *str1) {
        if (cur_depth == depth) {
            return;
        }

        T next_seeds[2];
        DLenPRG(seed, next_seeds);
        str0[cur_depth] ^= next_seeds[0];
        str1[cur_depth] ^= next_seeds[1];
        PrepareCorrelation(cur_depth + 1, depth, next_seeds[0], str0, str1);
        PrepareCorrelation(cur_depth + 1, depth, next_seeds[1], str0, str1);
    }

    void send_ot(const T* data0, const T* data1, int64_t length) {
		// block * data = new block[length];
        AlignedUnVector<std::array<block, 2>> data(length);
        sync_wait(sender->send(data, prg, chl));
        T *pad1 = new T[2*length];

        for (int64_t i = 0; i < length; ++i) {
            pad1[2*i] = *(T*)(data[i][0].data()) ^ data0[i];
            pad1[2*i+1] = *(T*)(data[i][1].data()) ^ data1[i];
        }
        // macoro::sync_wait(chl.sync());
        // cout << "send size : " << 2 * length << endl;
        coproto::span<T> a(pad1, 2*length);
        // cout << a.size() << endl;
        sync_wait(chl.send(a));
        // delete[] pad;
        delete[] pad1;
		// delete[] data;
	}

	void recv_ot(T* datar, const BitVector &r, int64_t length) {
        vector<block> data(length);

        sync_wait(recver->receive(r, data, prg, chl));
        // block *data = new block[length+ot_bsize];
		// ferretcot->recv_cot(data, r, length);
		// block s;
		// io->flush();
		// io->recv_block(&s,1);
		// io->flush();
		// mitccrh.setS(s);
		// io->flush();

		T *res = new T[2*length];
		// block *pad = new block[length+ot_bsize];
        // block *pad = data;
        // cout << "ot_data" << endl;
		// for(int64_t i = 0; i < length; i+=ot_bsize) {
		// 	memcpy(pad+i, data+i, min(ot_bsize,length-i)*sizeof(block));
		// 	mitccrh.hash<ot_bsize, 1>(pad+i);
		// }
		// io->flush();
        // macoro::sync_wait(chl.sync());
        // cout << "recv size : " << 2 * length << endl;
        coproto::span<T> a(res, 2*length);
        // cout << a.size() << endl;
        sync_wait(chl.recv(a));
        // cout << 1 << endl;
        // io->recv_data(res, 2*sizeof(T)*length);
		// io->flush();
        ////#pragma omp parallel for
        for(int64_t i = 0; i < length; i+=8) {
			for(int64_t j = 0; j < 8 and j < length-i; ++j) {
				datar[i+j] = res[2*(i+j)+r[i+j]] ^ *(T*)(data[i+j].data());
                // cout << hex << (uint64_t)datar[i+j] << " " << r[i+j] << endl;
			}
		}
        delete[] res;
        // delete[] pad;
        // delete[] data;
	}
};
