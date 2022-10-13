##10/04 release feedback
* On 04.10.2022 17:36, Lyu, Xuanyu wrote:
> Hi CRAN reviewers,
>
> I think those two notes you detected are both false positives.
>
>     1. Visscher is the last name of the author of the cited paper. It is not
>     a misspelled word.
>     2. The Description indeed contains complete sentences.

The last one is not a sentence, Write "See ..... ." (i. e. also end with
a period).
Fixed.

Please fix and resubmit.

## 09/30 release feedback

* Please always explain all acronyms in the description text. -> 'ACE'
Explained.

* Please do not start the description with "This is a package", package name, title or similar.
Done!

* If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
** authors (year) <doi:...>
** authors (year) <arXiv:...>
** authors (year, ISBN:...)
** or if those are not available: <https:...> 
**with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")
Added the citation. 

* In ypur LICENSE-file you claim "SimFit authors" to be the copyrightholders. Who are they? Or do you mean yourself? In that case please write "ACEsimFit authors" for better clarity or instead write their names.
Changed. 

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
