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

				  С++98 БИБЛИОТЕКА ДЛЯ РАБОТЫ С DE-ЭФЕМЕРИДАМИ
				   ВЕРСИЯ 0.4 [DEVELOP/UNPUBLISHED/C++BUILDER]

*******************************************************************************/

#include <fstream>
#include <cstring>
#include <stdint.h>
#include <string>
#include <map>
#include <vector>

/* Формат DE-эфемерид. */
#define DPH_RLS_LABELS_COUNT 3  // Число строк в метке выпуска.
#define DPH_RLS_LABEL_SIZE 84   // Длина строки метки.
#define DPH_CNAME_SIZE 6        // Длина имени константы.
#define DPH_CCOUNT_MAX_OLD 400  // Кол-во констант (стар.)
#define DPH_CCOUNT_MAX_NEW 1000 // Кол-во констант (нов.).

namespace dph {

/* -----------------------------| HELP VALUES |------------------------------ */

/* Значения параметров для метода dph::EphemerisReelase::calculateBody(...).
   Данные значения соответствуют индексам тел, для которых можно получить
   результат вычислений. */
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
	B_SSBARY  = 12, /* Барицентр Солнечной Системы. */
	B_EMBARY  = 13  /* Барицентр системы Земля-Луна. */
};

/* Значения параметров для метода dph::EphemerisRelease::calculateOther(...).
   Данные значения соответствуют индексам прочих элементов, для которых можно
   получить результат вычислений. */
enum Other
{
	O_EARTH_NUTATIONS               = 14,
	O_LUNAR_MANTLE_LIBRATION        = 15,
	O_LUNAR_MANTLE_ANGULAR_VELOCITY = 16,
	O_TTmTDB                        = 17
};

/* Значения параметров для методов dph::EphemerisRelease::calculateBody(...) и
   dph::EphemerisRelease::calculateOther(...).
   Данные значения соответствуют индексам результатов вычислений, которые можно
   получить. */
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

    /* [СТАНДАРТНЫЕ МЕТОДЫ] */

    /* Конструктор по умолчанию. */
    EphemerisRelease();

    /* Конструктор по пути к бинарному файлу эфемерид.
            Чтение файла, проверка полученных значений. */
    explicit EphemerisRelease(const std::string& binaryFilePath);

    /* Конструктор копирования.
            Проверка полученных значений и доступ к файлу.
            При неудачной проверке объект очищается. */
    EphemerisRelease(const EphemerisRelease& other);

    /* Оператор копирования.
            Проверка полученных значений и доступ к файлу.
            При неудачной проверке объект очищается. */
    EphemerisRelease& operator=(const EphemerisRelease& other);

    /* Деструктор.
            Просто деструктор.*/
    ~EphemerisRelease();


    /* [ПОЛЬЗОВАТЕЛЬСКИЕ МЕТОДЫ] */

    /* Открыть бинарый файл эфемерид. */
    void open(const std::string& binaryFilePath);

    /* Радиус-вектор (или вектор состояния) одного тела относительно другого. */
    void calculateBody(unsigned calculationResult, unsigned targetBody,
        unsigned centerBody, double JED, double* resultArray) const;

	/* Радиус-вектор (или вектор состояния) одного тела относительно другого. */
	void calculate(unsigned targetBody, unsigned centerBody,
		double JED, double r[3], double v[3]) const;

	/* Значение одного из дополнительных элементов, хранящихся в выпуске DE. */
    void calculateOther(unsigned calculationResult, unsigned otherItem,
        double JED, double* resultArray) const;

    /* Готовность объекта к использованию. */
    bool isReady() const;

    /* Первая доступная дата для рассчётов. */
    double startDate() const;

    /* Последняя доступная дата для рассчётов. */
    double endDate() const;

    /* Номер выпуска. */
    uint32_t releaseIndex() const;

    /* Заголовочная информация DE. */
    const std::string& releaseLabel() const;

    /* Значение константы по её имени. */
    double constant(const std::string& constantName) const;


private:
    /* Заголовочная информация выпуска DE. */
	std::string m_releaseLabel;  // Метка выпуска.
	uint32_t    m_releaseIndex;  // Номерная часть индекса выпуска.
	double      m_startDate;     // Дата начала выпуска (JED).
	double      m_endDate;       // Дата окончания выпуска (JED).
	double      m_blockTimeSpan; // Временная протяжённость блока.
	uint32_t    m_keys[15][3];   // Ключи поиска коэффициентов.
	double      m_au;            // Астрономическая единица [км].
	double      m_emrat;         // Отношение массы Земли к массе Луны.
	std::map<std::string, double> m_constants; // Константы выпуска.

    /* Производные значения для работы с выпуском DE. */
	size_t   m_blocksCount;     // Количество блоков в файле.
	uint32_t m_ncoeff;          // Количество коэффициентов в блоке.
	double   m_emrat2;          // Отношение массы Луны к массе Земля-Луна.
	double   m_dimensionFit;    // Значение для соблюдения размерности.
	size_t   m_blockSize_bytes; // Размер блока в байтах.

	/* Обеспечение работоспособности объекта. */
    bool m_ready;                             // Готовность объекта к работе.
    std::string m_binaryFilePath;             // Путь к файлу эфемерид.
    mutable std::ifstream m_binaryFileStream; // Поток чтения файла.
    mutable std::vector<double> m_buffer;	  // Буффер блока с коэффициентами.
    mutable std::vector<double> m_poly;		  // Значения полиномов.
    mutable std::vector<double> m_dpoly;	  // Значения производных полиномов.


    /* [ПРИВАТНЫЕ МЕТОДЫ] */

    /* Обрезать повторяющиеся пробелы (' ') с конца массива символов "charArray"
       размера "arraySize". */
    static std::string cutBackSpaces(const char* charArray, size_t arraySize);

    /* Приведение объекта к изначальному состоянию. */
    void clear();

    /* Копирование информации из объекта "other" в текущий объект. */
    void copyHere(const EphemerisRelease& other);


    /* [Работа с файлом] */

    /* Считать данные из бинарного файла выпуска DE в объект. */
    void readAndPackData();

    /* Дополнительные вычисления после чтения файла.
       Входит в состав readAndPackData(). */
    void additionalCalculations();

    /* Проверка значений, хранящихся в объекте и проверка файла. */
    bool isDataCorrect() const;

    /* Проверка начальных и конечных дат всех блоков в файле.
       Подтверждает целостность файла и доступность всех коэффициентов.
       Входит в состав проверки isDataCorrect(). */
    bool check_blocksDates() const;

    /* Заполнение буффера "m_buffer" коэффициентами требуемого блока. */
    void fillBuffer(size_t block_num) const;


    /* [Вычисления] */

    /* Интерполяция компонент базового элемента. */
    void interpolatePosition(unsigned baseItemIndex, double normalizedTime,
        const double* coeffArray, unsigned componentsCount,
        double* resultArray) const;

    /* Интерполяция компонент и их первых производных базового элемента. */
    void interpolateState(unsigned baseItemIndex, double normalizedTime,
        const double* coeffArray, unsigned componentsCount,
        double* resultArray) const;

    /* Значения компонент базового элемента. */
    void calculateBaseItem(unsigned baseItemIndex, double JED,
        unsigned calculationResult , double* resultArray) const;

    /* Радиус-вектор (или вектор состояния) Земли относительно барицентра СС. */
    void calculateBaseEarth(double JED, unsigned calculationResult,
        double* resultArray) const;

    /* Радиус-вектор (или вектор состояния) Луны относительно барицентра СС. */
    void calculateBaseMoon(double JED, unsigned calculationResult,
        double* resultArray) const;

}; /* class EphemerisRelease */
/* -------------------------------------------------------------------------- */

} /* namespace dph */

//---------------------------------------------------------------------------
#endif
