/**
 * @file CPolynomialGenerator.h
 * @brief Class definition for polynomials generator
 * @author Oren Cohen
 * @bug No known bugs
 */

#ifndef LOOKUP_TABLES_CPOLYNOMIALGENERATOR_H_
#define LOOKUP_TABLES_CPOLYNOMIALGENERATOR_H_

#include <string>

#include <list>
#include <vector>

/**
 * Generates polynomials
 *
 * Usage example:
 * @code
 * CPolynomialGenerator generator(2, 3);
 * This will generates : 101 111
 * @endcode
 */
class CPolynomialGenerator
{
public:

	/**
	 * Parameterized constructor
	 * @param p_nN Number of polynomials needed
	 * @param p_nK Polynomial symbol size, meaning the degree of the polynomial
	 */
	CPolynomialGenerator(unsigned int const p_cnN, unsigned int const p_cnK);

	/**
	 * Destructor
	 */
	virtual ~CPolynomialGenerator(void);

	/**
	 * Get the degree of the polynomial generator
	 * @return The degree of the polynomials
	 */
	unsigned int GetDegree(void) const;
	/**
	 * Get the number of polynomials that generated
	 * @return The number of polynomials
	 */
	unsigned int GetSize(void) const;
	/**
	 * Get specific polynomial
	 * @param p_nIndex The index of the polynomial
	 * @return Chosen polynomial
	 */
	unsigned int GetPolynomial(unsigned int const p_cnIndex) const;

	/**
	 * For Debug purposes - Generate binary string for an integer
	 * @param p_nNumber The converted number
	 * @param p_nLength The binary string desired size
	 * @return Binary string at specific size
	 */
	static std::string const GetBinaryString(unsigned int const p_cnNumber, unsigned int const p_cnLength);
	static std::string const GetBinaryString(std::list<unsigned int> const& p_rccNumber, unsigned int const p_cnLength);


	/**
	 * For Debug purposes - overload the << operator for this class
	 * @param p_rcOs The output stream
	 * @param p_rcGenerator The polynomials generator for output
	 * @return The current output stream
	 */
	friend std::ostream& operator<<(std::ostream& p_rcOs, CPolynomialGenerator const& p_rccGenerator);



private:

	/**
	 * Create polynomials
	 */
	void CreatePolynomials(void);

	/**
	 * The degree of the polynomials
	 */
	unsigned int const m_cnDegree;
	/**
	 * Polynomials buffer
	 */
	std::vector<unsigned int> m_cPolynomials;

};

#endif /* LOOKUP_TABLES_CPOLYNOMIALGENERATOR_H_ */
