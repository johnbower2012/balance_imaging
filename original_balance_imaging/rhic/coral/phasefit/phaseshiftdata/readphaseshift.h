//ReadPhaseShift.h
#ifndef	READPHASESHIFT_H
#define	READPHASESHIFT_H

class	CReadPhaseShift
{
public:
	char szInteraction[4];
	char szPartialWave[3];
	double	ReadPhaseShift(double	q);
	CReadPhaseShift(const char szInteraction[], const char szPartialWave[]);

protected:
	double	m1, m2;	//m2 is defined as the moving particle and m1 is the target
	bool	ReadData();
	const static int	MAX_DATA = 50;
	double	plab[MAX_DATA];
	double	del[MAX_DATA];

};

#endif	//READPHASESHIFT_H