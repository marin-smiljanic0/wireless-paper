/*
 * Input:
 * n - number of subcarriers
 * k - number of chunks
 * P_sr - total power on SR
 * P_rd - total power on RD
 * pi - n numbers, power allocation on SR
 * gi - n numbers, channel gain on SR
 * hi - n numbers, channel gain on RD
 * ai - k numbers, ratios of data rates
 */

#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <cassert>
#include <cmath>

double q_function(double x) {
    return 0.5 * erfc(x / sqrt(2));
}

double calc_ber(int r, double gamma) {
    if (r == 0)
	return 0.0;
    return 4.0 / r * (1.0 - pow(2.0, -r / 2.0)) * q_function(sqrt(3.0 * gamma / (pow(2.0, r) - 1.0)));
}

double square(double x) {
    return x;
    //    return x * x;
}

double log2(double x) {
    return log(x) / log(2);
}

class subchannel {

public:

    double getSnr() const {
	return snr_;
    }

    double getPower() const {
	return power_;
    }

    double getH() const {
	return h_;
    }

    int getIndex() const {
	return index_;
    }

    subchannel(double power, double h, int index) :
	power_(power),
	h_(h),
	index_(index),
	snr_(power * square(h)) {}

    subchannel(const subchannel &other) :
	power_(other.power_),
	h_(other.h_),
	index_(other.index_),
	snr_(other.snr_) {}

    bool operator < (const subchannel &other) const {
	if (snr_ != other.snr_) {
	    return snr_ > other.snr_;
	}
	return index_ < other.index_;
    }

    /*
     * Joint SNR of a subsystem.
     * snr1 - the SNR of the SR subchannel
     * snr2 - the SNR of the RD subchannel
     */
    double match_snr(const subchannel &other) const {
	return snr_ * other.snr_ / (snr_ + other.snr_ + 1);
    }

private:
    double power_;
    double h_;
    int index_;
    double snr_;
    
};

struct chunk {

public:    
    double getAverageSnr() const {
	return average_snr_;
    }

    std::pair<int, int> getRange() const {
	return std::pair<int, int>(subchannels_[0].getIndex(), subchannels_.back().getIndex());
    }

    std::vector<subchannel> getSubchannels() const {
	return subchannels_;
    }

    // [left, right>
    chunk(const std::vector<subchannel> &v, int left, int right) {
	subchannels_.reserve(right - left);

	double sum = 0;
	for (int i = left; i < right; ++i) {
	    subchannels_.push_back(v[i]);
	    sum += v[i].getSnr();
	}
	average_snr_ = sum / (right - left);
    }
    
    bool operator < (const chunk &other) const {
	if (average_snr_ != other.average_snr_) {
	    return average_snr_ > other.average_snr_;
	}
	return subchannels_[0] < other.subchannels_[0];
    }

private:
    std::vector<subchannel> subchannels_;
    double average_snr_;
};

namespace Input {
    template<typename T> int scanf_t(FILE *fin, T *addr) { return 0; }

    template<>
    int scanf_t<int>(FILE *fin, int *addr) {
	return fscanf(fin, "%d", addr);
    }

    template<>
    int scanf_t<double>(FILE *fin, double *addr) {
	return fscanf(fin, "%lf", addr);
    }

    template<typename T> std::vector<T> read_vector(FILE *fin, size_t n) {
	std::vector<T> v;
	v.reserve(n);

	for (size_t i = 0; i < n; ++i) {
	    T x;
	    scanf_t<T>(fin, &x);
	    v.push_back(x);
	}

	return v;
    }
};

namespace Output {
    template<typename T> int printf_t(T data) { return 0; }

    template<>
    int printf_t<int>(int data) {
	return printf("%d ", data);
    }

    template<>
    int printf_t<double>(double data) {
	return printf("%.6lf ", data);
    }

    template<>
    int printf_t<subchannel>(subchannel data) {
	return printf("%d %.2lf ", data.getIndex(), data.getSnr());
    }

    template<>
    int printf_t<chunk>(chunk data) {
    	return printf("%d %d %.2lf ",
    		      data.getRange().first,
    		      data.getRange().second,
    		      data.getAverageSnr());
    }

    template<typename T> void print_vector(const std::vector<T> &v) {
	size_t n = v.size();
	for (size_t i = 0; i < n; ++i) {
	    printf_t(v[i]);
	}
	putchar('\n');
    }
};

std::vector<chunk> make_chunks(const std::vector<subchannel> &subchannels, int k) {
    std::vector<chunk> ret;
    ret.reserve(k);
    int chunk_size = subchannels.size() / k;

    for (size_t pos = 0; pos < subchannels.size(); pos += chunk_size) {
	ret.push_back(chunk(subchannels, pos, pos + chunk_size));
    }
    return ret;
}

std::vector<double> get_powers(const std::vector<subchannel> &sr_subchannels,
			       const std::vector<subchannel> &rd_subchannels,
			       const std::vector<double> &g,
			       const std::vector<double> &h,
			       const std::vector<double> &a,
			       double p_total,
			       double l) {
    size_t n = sr_subchannels.size();
    size_t n_chunks = a.size();
    size_t chunk_size = n / n_chunks;

    std::vector<double> ret(n);

    for (size_t i = 0; i < n; ++i) {
	double ratio = a[i / chunk_size];
	double num = (pow(2.0, l * ratio / chunk_size) - 1)
	    * (sr_subchannels[i].getPower() * square(sr_subchannels[i].getH()) + 1);
	double denom = square(rd_subchannels[i].getH())
	    * (sr_subchannels[i].getPower() * square(sr_subchannels[i].getH())
	       - pow(2.0, l * ratio / chunk_size) + 1);
	ret[rd_subchannels[i].getIndex()] = num / denom;
    }
    
    return ret;
}

double achievable(const std::vector<subchannel> &sr_subchannels,
		  const std::vector<subchannel> &rd_subchannels,
		  const std::vector<double> &g,
		  const std::vector<double> &h,
		  const std::vector<double> &a,
		  double p_total,
		  double l) {
    size_t n = sr_subchannels.size();
    size_t n_chunks = a.size();
    size_t chunk_size = n / n_chunks;

    double ret = 0;

    for (size_t i = 0; i < n; ++i) {
	double ratio = a[i / chunk_size];
	double num = (pow(2.0, l * ratio / chunk_size) - 1)
	    * (sr_subchannels[i].getPower() * square(sr_subchannels[i].getH()) + 1);
	double denom = square(rd_subchannels[i].getH())
	    * (sr_subchannels[i].getPower() * square(sr_subchannels[i].getH())
	       - pow(2.0, l * ratio / chunk_size) + 1);
	ret += num / denom;
    }

    return ret;
}

std::vector<subchannel> get_subchannels(const std::vector<chunk> &v) {
    std::vector<subchannel> ret;
    ret.reserve(v.size() * v[0].getSubchannels().size());

    for (size_t i = 0; i < v.size(); ++i) {
	for (size_t j = 0; j < v[i].getSubchannels().size(); ++j) {
	    ret.push_back(v[i].getSubchannels()[j]);
	}
    }

    return ret;
}

void print_vector_chunk(const std::vector<chunk> &v) {
    for (size_t i = 0; i < v.size(); ++i) {
	printf("%d %d %.2lf ... ", v[i].getRange().first, v[i].getRange().second,
	       v[i].getAverageSnr());
	for (size_t j = 0; j < v[i].getSubchannels().size(); ++j) {
	    printf("%.2lf ", v[i].getSubchannels()[j].getSnr());
	}
	putchar('\n');
    }
}

double data_rate(double snr) {
    return std::max(0.0, log2(1 + snr));
}

std::vector<double> get_ratios(int chunk_size, const std::vector<subchannel> &sr, const std::vector<subchannel> &rd, const std::vector<double> &userRatios) {
    std::vector<double> chunkRates;
    double sumRates = 0.0;
    int n = (int) sr.size();

    for (int i = 0; i < n / chunk_size; ++i) {
	double rate = 0.0;
	for (int j = 0; j < chunk_size; ++j) {
	    rate += data_rate(sr[i * chunk_size + j].match_snr(rd[i * chunk_size + j]));
	}
	chunkRates.push_back(rate);
	sumRates += rate;
    }

    double sumUsers = 0.0;
    for (int i = 0; i < userRatios.size(); ++i) {
	sumUsers += userRatios[i];
    }

    std::vector<double> req;
    std::vector<double> asgn(userRatios.size());
    fill(asgn.begin(), asgn.end(), 0.0);

    for (int i = 0; i < userRatios.size(); ++i) {
	req.push_back(userRatios[i] / sumUsers * sumRates);
    }

    std::vector<int> ind(n / chunk_size);
    fill(ind.begin(), ind.end(), -1);
    for (int i = 0; i < chunkRates.size(); ++i) {
	double max_diff = 0, max_ind = -1;
	for (int j = 0; j < req.size(); ++j) {
	    if (req[j] - asgn[j] > max_diff) {
		max_diff = req[j] - asgn[j];
		max_ind = j;
	    }
	}
	asgn[max_ind] += chunkRates[i];
	ind[i] = max_ind;
    }

    std::vector<double> final_ratios;
    for (int i = 0; i < chunkRates.size(); ++i) {
	final_ratios.push_back(chunkRates[i] / asgn[ind[i]] * userRatios[ind[i]]);
    }
    return final_ratios;
}

std::pair<double, double> run(int n, int n_chunks, double p_sr, double p_rd, const std::vector<double> &ratios, FILE *ch1, FILE *ch2) {
    std::vector<double> g = Input::read_vector<double>(ch1, n);
    std::vector<double> h = Input::read_vector<double>(ch2, n);

    // for (int i = 0; i < n; ++i)
    // 	printf("%.2lf %.2lf\n", g[i], h[i]);

    std::vector<subchannel> snr_sr;
    snr_sr.reserve(n);
    std::vector<subchannel> snr_rd;
    snr_rd.reserve(n);
    for (int i = 0; i < n; ++i) {
	snr_sr.push_back(subchannel(p_sr / n, g[i], i));
	snr_rd.push_back(subchannel(p_rd / n, h[i], i));
    }

    // for (int i = 0; i < n; ++i)
    // 	printf("%d\n", data_rate(snr_sr[i].match_snr(snr_rd[i])));
    
    std::vector<chunk> sr_chunks = make_chunks(snr_sr, n_chunks);
    std::vector<chunk> rd_chunks = make_chunks(snr_rd, n_chunks);
    
    sort(sr_chunks.begin(), sr_chunks.end());
    sort(rd_chunks.begin(), rd_chunks.end());

    for (int i = 0; i < sr_chunks.size(); ++i) {
	std::pair<int, int> r1 = sr_chunks[i].getRange();
	std::pair<int, int> r2 = rd_chunks[i].getRange();
	// printf("%d %d -----> %d %d\n", r1.first, r1.second, r2.first, r2.second);
	// printf("%.2lf %.2lf\n", sr_chunks[i].getAverageSnr(), rd_chunks[i].getAverageSnr());
    }

    std::vector<subchannel> sr_subchannels = get_subchannels(sr_chunks);
    std::vector<subchannel> rd_subchannels = get_subchannels(rd_chunks);

    std::vector<double> a_int = get_ratios(n / n_chunks, sr_subchannels, rd_subchannels, ratios);
    std::vector<double> a;
    double sum = 0;
    for (int i = 0; i < n_chunks; ++i) {
    	sum += a_int[i];
    }
    for (int i = 0; i < n_chunks; ++i) {
    	a.push_back(1.0 * a_int[i] / sum);
    }
    sort(a.begin(), a.end(), std::greater<double>());

    // for (int i = 0; i < n; ++i)
    // 	printf("%d\n", data_rate(sr_subchannels[i].match_snr(rd_subchannels[i])));
    // printf("-----------------\n");

    double max_lambda = 10000000.0;
    size_t worst = -1;
    for (size_t i = 0; i < sr_subchannels.size(); ++i) {
	double ratio = a[i / (n / n_chunks)];
	double denom = sr_subchannels[i].getPower() * square(sr_subchannels[i].getH()) + 1;

	// 2 ^ (ratio * lambda / chunk_size) = denom
	double lambda = log2(denom) * (n / n_chunks) / ratio;
	max_lambda = std::min(lambda, max_lambda);
	if (fabs(lambda - max_lambda) < 1e-9) {
	    worst = i;
	}
    }

    double lo = 0, hi = max_lambda;
    while (hi - lo > 1e-6) {
    	double mid = (lo + hi) / 2.0;
    	double ach = achievable(sr_subchannels, rd_subchannels, g, h, a, p_rd, mid);
    	if (ach < p_rd) {
    	    lo = mid;
    	} else {
    	    hi = mid;
    	}
    }

    std::vector<double> powers = get_powers(sr_subchannels, rd_subchannels, g, h, a, p_rd, lo);
    double p_sum = 0.0;
    for (size_t i = 0; i < powers.size(); ++i) {
    	p_sum += powers[i];
    }

    double rate_sum = 0.0;
    double ber_sum = 0.0;
    int chunk_size = n / n_chunks;
    for (int i = 0; i < n; ++i) {
	if ((int) i % chunk_size == 0) {
	    std::pair<int, int> range_sr = sr_chunks[i / chunk_size].getRange();
	    std::pair<int, int> range_rd = rd_chunks[i / chunk_size].getRange();
	}
	subchannel rd_new(powers[rd_subchannels[i].getIndex()], rd_subchannels[i].getH(), i);
	double curr_snr = rd_new.match_snr(sr_subchannels[i]);
	double curr_rate = data_rate(curr_snr);
	rate_sum += curr_rate;
	ber_sum = calc_ber((int) curr_rate, curr_snr);
    }

    return std::pair<double, double>(rate_sum, ber_sum / n);
}


int main(int argc, char **argv) {
    int n = atoi(argv[1]);
    double p_sr = atof(argv[2]);
    double p_rd = atof(argv[3]);

    std::vector<double> userRatios;
    userRatios.push_back(1);
    userRatios.push_back(2);
    userRatios.push_back(4);
    userRatios.push_back(8);
    std::cout << "FINAL\n";

    for (int n_relays = 1; n_relays <= 1; ++n_relays) {
	for (int chunks = 4; chunks <= n; chunks *= 2) {
	    double sum = 0;
	    double ber_sum = 0.0;
	    //for (int i = 1; i <= 100; ++i) {
	    //for (int j = 1; j <= 100; ++j) {
	    for (int i = 1; i <= 20; ++i) {
		for (int j = 1; j <= 20; ++j) {
		    //std::cout << chunks << " " << i*100 + j << std::endl;
		    std::ostringstream s;
		    s << "test_data_new_large/kanal" << i << ".txt";
		    std::string s1 = s.str();
		    FILE *fin = fopen(s1.c_str(), "r");

		    std::ostringstream s2;
		    s2 << "test_data_new_large/kanal" << j << ".txt";
		    std::string s3 = s2.str();
		    FILE *fin2 = fopen(s3.c_str(), "r");
		    std::pair<double, double> x = run(n, chunks, p_sr, p_rd, userRatios, fin, fin2);
		    sum += x.first;
		    ber_sum += x.second;
		    fclose(fin);
		    fclose(fin2);
		}
	    }
	    printf("(%d %.2lf) ", chunks, 1.0 * sum / 400);
	}

	printf("\n");
    }

    return 0;
}
