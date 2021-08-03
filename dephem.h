//---------------------------------------------------------------------------

#ifndef dephemH
#define dephemH

/*******************************************************************************

      8888888b.  8888888888 8888888b.  888    888 8888888888 888b     d888
      888  "Y88b 888        888   Y88b 888    888 888        8888b   d8888
      888    888 888        888    888 888    888 888        88888b.d88888
      888    888 8888888    888   d88P 8888888888 8888888    888Y88888P888
      888    888 888        8888888P"  888    888 888        888 Y888P 888
      888    888 888        888        888    888 888        888  Y8P  888
      888  .d88P 888        888        888    888 888        888   "   888
      8888888P"  8888888888 888        888    888 8888888888 888       888

				  �++98 ���������� ��� ������ � DE-�����������
				   ������ 0.4 [DEVELOP/UNPUBLISHED/C++BUILDER]

*******************************************************************************/

#include <fstream>
#include <cstring>
#include <stdint.h>
#include <string>
#include <map>
#include <vector>

/* ������ DE-��������. */
#define DPH_RLS_LABELS_COUNT 3  // ����� ����� � ����� �������.
#define DPH_RLS_LABEL_SIZE 84   // ����� ������ �����.
#define DPH_CNAME_SIZE 6        // ����� ����� ���������.
#define DPH_CCOUNT_MAX_OLD 400  // ���-�� �������� (����.)
#define DPH_CCOUNT_MAX_NEW 1000 // ���-�� �������� (���.).

namespace dph {

/* -----------------------------| HELP VALUES |------------------------------ */

/* �������� ���������� ��� ������ dph::EphemerisReelase::calculateBody(...).
   ������ �������� ������������� �������� ���, ��� ������� ����� ��������
   ��������� ����������. */
enum Body
{
	B_MERCURY = 1,
	B_VENUS   = 2,
	B_EARTH   = 3,
	B_MARS    = 4,
	B_JUPITER = 5,
	B_SATURN  = 6,
	B_URANUS  = 7,
	B_NEPTUNE = 8,
	B_PLUTO   = 9,
	B_MOON    = 10,
	B_SUN     = 11,
	B_SSBARY  = 12, /* ��������� ��������� �������. */
	B_EMBARY  = 13  /* ��������� ������� �����-����. */
};

/* �������� ���������� ��� ������ dph::EphemerisRelease::calculateOther(...).
   ������ �������� ������������� �������� ������ ���������, ��� ������� �����
   �������� ��������� ����������. */
enum Other
{
	O_EARTH_NUTATIONS               = 14,
	O_LUNAR_MANTLE_LIBRATION        = 15,
	O_LUNAR_MANTLE_ANGULAR_VELOCITY = 16,
	O_TTmTDB                        = 17
};

/* �������� ���������� ��� ������� dph::EphemerisRelease::calculateBody(...) �
   dph::EphemerisRelease::calculateOther(...).
   ������ �������� ������������� �������� ����������� ����������, ������� �����
   ��������. */
enum Calc
{
	CALC_POS = 0,
	CALC_STATE = 1
};

/* -------------------------------------------------------------------------- */


/* -----------------------| CLASS EPHEMERIS_RELEASE |------------------------ */
class EphemerisRelease
{
public:

    /* [����������� ������] */

    /* ����������� �� ���������. */
    EphemerisRelease();

    /* ����������� �� ���� � ��������� ����� ��������.
            ������ �����, �������� ���������� ��������. */
    explicit EphemerisRelease(const std::string& binaryFilePath);

    /* ����������� �����������.
            �������� ���������� �������� � ������ � �����.
            ��� ��������� �������� ������ ���������. */
    EphemerisRelease(const EphemerisRelease& other);

    /* �������� �����������.
            �������� ���������� �������� � ������ � �����.
            ��� ��������� �������� ������ ���������. */
    EphemerisRelease& operator=(const EphemerisRelease& other);

    /* ����������.
            ������ ����������.*/
    ~EphemerisRelease();


    /* [���������������� ������] */

    /* ������� ������� ���� ��������. */
    void open(const std::string& binaryFilePath);

    /* ������-������ (��� ������ ���������) ������ ���� ������������ �������. */
    void calculateBody(unsigned calculationResult, unsigned targetBody,
        unsigned centerBody, double JED, double* resultArray) const;

	/* ������-������ (��� ������ ���������) ������ ���� ������������ �������. */
	void calculate(unsigned targetBody, unsigned centerBody,
		double JED, double r[3], double v[3]) const;

	/* �������� ������ �� �������������� ���������, ���������� � ������� DE. */
    void calculateOther(unsigned calculationResult, unsigned otherItem,
        double JED, double* resultArray) const;

    /* ���������� ������� � �������������. */
    bool isReady() const;

    /* ������ ��������� ���� ��� ���������. */
    double startDate() const;

    /* ��������� ��������� ���� ��� ���������. */
    double endDate() const;

    /* ����� �������. */
    uint32_t releaseIndex() const;

    /* ������������ ���������� DE. */
    const std::string& releaseLabel() const;

    /* �������� ��������� �� � �����. */
    double constant(const std::string& constantName) const;


private:
    /* ������������ ���������� ������� DE. */
	std::string m_releaseLabel;  // ����� �������.
	uint32_t    m_releaseIndex;  // �������� ����� ������� �������.
	double      m_startDate;     // ���� ������ ������� (JED).
	double      m_endDate;       // ���� ��������� ������� (JED).
	double      m_blockTimeSpan; // ��������� ������������ �����.
	uint32_t    m_keys[15][3];   // ����� ������ �������������.
	double      m_au;            // ��������������� ������� [��].
	double      m_emrat;         // ��������� ����� ����� � ����� ����.
	std::map<std::string, double> m_constants; // ��������� �������.

    /* ����������� �������� ��� ������ � �������� DE. */
	size_t   m_blocksCount;     // ���������� ������ � �����.
	uint32_t m_ncoeff;          // ���������� ������������� � �����.
	double   m_emrat2;          // ��������� ����� ���� � ����� �����-����.
	double   m_dimensionFit;    // �������� ��� ���������� �����������.
	size_t   m_blockSize_bytes; // ������ ����� � ������.

	/* ����������� ����������������� �������. */
    bool m_ready;                             // ���������� ������� � ������.
    std::string m_binaryFilePath;             // ���� � ����� ��������.
    mutable std::ifstream m_binaryFileStream; // ����� ������ �����.
    mutable std::vector<double> m_buffer;	  // ������ ����� � ��������������.
    mutable std::vector<double> m_poly;		  // �������� ���������.
    mutable std::vector<double> m_dpoly;	  // �������� ����������� ���������.


    /* [��������� ������] */

    /* �������� ������������� ������� (' ') � ����� ������� �������� "charArray"
       ������� "arraySize". */
    static std::string cutBackSpaces(const char* charArray, size_t arraySize);

    /* ���������� ������� � ������������ ���������. */
    void clear();

    /* ����������� ���������� �� ������� "other" � ������� ������. */
    void copyHere(const EphemerisRelease& other);


    /* [������ � ������] */

    /* ������� ������ �� ��������� ����� ������� DE � ������. */
    void readAndPackData();

    /* �������������� ���������� ����� ������ �����.
       ������ � ������ readAndPackData(). */
    void additionalCalculations();

    /* �������� ��������, ���������� � ������� � �������� �����. */
    bool isDataCorrect() const;

    /* �������� ��������� � �������� ��� ���� ������ � �����.
       ������������ ����������� ����� � ����������� ���� �������������.
       ������ � ������ �������� isDataCorrect(). */
    bool check_blocksDates() const;

    /* ���������� ������� "m_buffer" �������������� ���������� �����. */
    void fillBuffer(size_t block_num) const;


    /* [����������] */

    /* ������������ ��������� �������� ��������. */
    void interpolatePosition(unsigned baseItemIndex, double normalizedTime,
        const double* coeffArray, unsigned componentsCount,
        double* resultArray) const;

    /* ������������ ��������� � �� ������ ����������� �������� ��������. */
    void interpolateState(unsigned baseItemIndex, double normalizedTime,
        const double* coeffArray, unsigned componentsCount,
        double* resultArray) const;

    /* �������� ��������� �������� ��������. */
    void calculateBaseItem(unsigned baseItemIndex, double JED,
        unsigned calculationResult , double* resultArray) const;

    /* ������-������ (��� ������ ���������) ����� ������������ ���������� ��. */
    void calculateBaseEarth(double JED, unsigned calculationResult,
        double* resultArray) const;

    /* ������-������ (��� ������ ���������) ���� ������������ ���������� ��. */
    void calculateBaseMoon(double JED, unsigned calculationResult,
        double* resultArray) const;

}; /* class EphemerisRelease */
/* -------------------------------------------------------------------------- */

} /* namespace dph */

//---------------------------------------------------------------------------
#endif
