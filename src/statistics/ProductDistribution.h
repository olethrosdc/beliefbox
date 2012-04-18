// -*- Mode: c++ -*-
#include <vector>

template <class P, class T>
class ProductDistribution
{
protected:
    std::vector<P*> distributions;
public:
    ProductDistribution(std::vector<P*> distributions_)
        : distributions(distributions_)
    {
    }

    virtual ~ProductDistribution() {}

    virtual std::vector<T> generate() const
    {
        std::vector<T> x(distributions.size());
        for (uint i=0; i<distributions.size(); ++i) {
            x[i] = distributions[i]->generate();
        }
        return x;
    }

    virtual real log_pdf(const std::vector<T>& x) const
    {
        real log_p = 0;
        for (uint i=0; i < distributions.size(); ++i) {
            log_p += distributions[i]->log_pdf(x[i]);
        }
        return log_p;
    }

    virtual real pdf(const std::vector<T>& x) const
    {
        return exp(log_pdf(x));
    }
};
