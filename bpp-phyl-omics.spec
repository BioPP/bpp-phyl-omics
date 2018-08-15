%define _prefix /usr

URL: https://github.com/BioPP/bpp-phyl-omics

Name: bpp-phyl-omics
Version: 2.4.1
Release: 1%{?dist}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: %{name}-%{version}.tar.gz
Summary: Bio++ Sequence library: genomics components
Group: Development/Libraries/C and C++
Requires: bpp-core = %{version}
Requires: bpp-seq = %{version}
Requires: bpp-phyl = %{version}
Requires: bpp-seq-omics = %{version}

BuildRoot: %{_builddir}/%{name}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: libbpp-core4 = %{version}
BuildRequires: libbpp-core-devel = %{version}
BuildRequires: libbpp-seq12 = %{version}
BuildRequires: libbpp-seq-devel = %{version}
BuildRequires: libbpp-phyl12 = %{version}
BuildRequires: libbpp-phyl-devel = %{version}
BuildRequires: libbpp-seq-omics3 = %{version}
BuildRequires: libbpp-seq-omics-devel = %{version}

AutoReq: yes
AutoProv: yes

%description
This library contains the genomics components of the Bio++ phylogenetics library.
It is part of the Bio++ project.

%package -n libbpp-phyl-omics3
Summary: Bio++ Phylogenetics library: genomics components
Group: Development/Libraries/C and C++

%description -n libbpp-phyl-omics3
This library contains the genomics components of the Bio++ phylogenetics library.
It is part of the Bio++ project.

%package -n libbpp-phyl-omics-devel
Summary: Bio++ Phylogenetics library: genomics components
Group: Development/Libraries/C and C++
Requires: libbpp-phyl-omics3 = %{version}
Requires: libbpp-seq12 = %{version}
Requires: libbpp-seq-devel = %{version}
Requires: libbpp-core4 = %{version}
Requires: libbpp-core-devel = %{version}
Requires: libbpp-phyl12 = %{version}
Requires: libbpp-phyl-devel = %{version}
Requires: libbpp-seq-omics3 = %{version}
Requires: libbpp-seq-omics-devel = %{version}

%description -n libbpp-phyl-omics-devel
The libbpp-phyl-omics-devel package contains the header files and static libraries for
building applications which use %{name}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DBUILD_TESTING=OFF"
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -n libbpp-phyl-omics3 -p /sbin/ldconfig

%postun -n libbpp-phyl-omics3 -p /sbin/ldconfig

%files -n libbpp-phyl-omics3
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so.*

%files -n libbpp-phyl-omics-devel
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%dir %{_prefix}/%{_lib}/cmake/
%dir %{_prefix}/%{_lib}/cmake/bpp-phyl-omics
%{_prefix}/%{_lib}/lib*.so
%{_prefix}/%{_lib}/lib*.a
%{_prefix}/%{_lib}/cmake/bpp-phyl-omics/bpp-phyl-omics*.cmake
%{_prefix}/include/*

%changelog
* Wed Aug 13 2018 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.4.1-1
- Compatibility update gcc8
* Mon Mar 12 2018 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.4.0-1
- Increased interface number
- Removed dynamic exceptions specifications.
* Tue Jun 06 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.1-1
- Increased interface number
* Wed May 10 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.0-1
- Added distance matrix output filter.
- Added clock and reparametrization options in model fitting.
* Fri Sep 26 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.2.0-1
- Added model fitting and parameter estimations.
* Thu Mar 07 2013 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.1.0-1
- Initial release. Contains tools to build phylogenies along a genome alignment (distance methods only for now).
* Tue Nov 06 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.3-1
- First draft of the spec file.

