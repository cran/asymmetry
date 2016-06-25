\name{studentmigration}
\alias{studentmigration}
\docType{data}
\title{
Student Migration in the Erasmus Program
}
\description{
The table lists the home and destination country of 268.142 students in the academic year 2012-2013. The 33 rows of this table refer to the home country whereas the 33 columns refer to the destination countries. The entries in the table are read as follows: 32 students from Bulgaria studied in The Netherlands, 18 students from the Netherlands studied in Bulgaria.
Macedonia (MK) was excluded from the published table because only one student from Macedonia studie abroad and this country did not receive any students.
}
\usage{data(studentmigration)}
\format{
 A matrix of 33 rows by 33 columns
}
\details{
The Erasmus program is a student exchange program from the European Union. Three million students had taken part since the start of the program in 1987. To join the program a student has study at least three months or do an internship of at least two months in another country.
The 2-letter codes shown below are supplied by the ISO (International Organization for Standardization). You can find them here: \href{http://publications.europa.eu/code/pdf/370000en.htm}{Countrycodes}
}
\source{
\href{http://ec.europa.eu/education/library/statistics/ay-12-13/annex-2_en.pdf}{http://ec.europa.eu/education/library/statistics/ay-12-13/annex-2_en.pdf}
}
\note{ Macedonia has been removed from this table because only one student from
this country participated in the program.}
\examples{
data(studentmigration)
hmap(studentmigration)
}
\keyword{datasets}
