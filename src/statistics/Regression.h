/** This is an abstract class for regression */
class Regression
{
public:
    /// Destructor
    virtual ~Regression()
    {
    }
    /// Get the expected value for a given x
    virtual Vector ExpectedValue(const Vector& x) = 0;
    /// Return the probability density at a specific y
    virtual real pdf(const Vectort& x, const Vector& y) = 0;
	/// Observe a particular input and output pair, returning the probability
    virtual real Observe(const Vectort& x, const Vector& y) = 0;
    /// Observe a particular input and class distribution
	/// Observe a particular input and output pair
    virtual real AddElement(const Vectort& x, const Vector& y) = 0;
};

