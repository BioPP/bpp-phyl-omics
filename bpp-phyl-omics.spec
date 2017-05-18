%define _basename bpp-phyl-omics
%define _version 2.3.0
%define _release 1
%define _prefix /usr

URL: http://biopp.univ-montp2.fr/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: http://biopp.univ-montp2.fr/repos/sources/%{_basename}-%{_version}.tar.gz
Summary: Bio++ Sequence library: genomics components
Group: Development/Libraries/C and C++
Requires: bpp-core = %{_version}
Requires: bpp-seq = %{_version}
Requires: bpp-phyl = %{_version}
Requires: bpp-seq-omics = %{_version}

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.8.11
BuildRequires: gcc-c++ >= 4.7.0
BuildRequires: libbpp-core2 = %{_version}
BuildRequires: libbpp-core-devel = %{_version}
BuildRequires: libbpp-seq9 = %{_version}
BuildRequires: libbpp-seq-devel = %{_version}
BuildRequires: libbpp-phyl9 = %{_version}
BuildRequires: libbpp-phyl-devel = %{_version}
BuildRequires: libbpp-seq-omics1 = %{_version}
BuildRequires: libbpp-seq-omics-devel = %{_version}

AutoReq: yes
AutoProv: yes

%description
This library contains the genomics components of the Bio++ phylogenetics library.
It is part of the Bio++ project.

%package -n libbpp-phyl-omics1
Summary: Bio++ Phylogenetics library: genomics components
Group: Development/Libraries/C and C++

%description -n libbpp-phyl-omics1
This library contains the genomics components of the Bio++ phylogenetics library.
It is part of the Bio++ project.

%package -n libbpp-phyl-omics-devel
Summary: Bio++ Phylogenetics library: genomics components
Group: Development/Libraries/C and C++
Requires: libbpp-phyl-omics1 = %{_version}
Requires: libbpp-seq9 = %{_version}
Requires: libbpp-seq-devel = %{_version}
Requires: libbpp-core2 = %{_version}
Requires: libbpp-core-devel = %{_version}
Requires: libbpp-phyl9 = %{_version}
Requires: libbpp-phyl-devel = %{_version}
Requires: libbpp-seq-omics1 = %{_version}
Requires: libbpp-seq-omics-devel = %{_version}

%description -n libbpp-phyl-omics-devel
The libbpp-phyl-omics-devel package contains the header files and static libraries for
building applications which use %{_basename}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DBUILD_TESTING=OFF"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -n libbpp-phyl-omics1 -p /sbin/ldconfig

%postun -n libbpp-phyl-omics1 -p /sbin/ldconfig

%files -n libbpp-phyl-omics1
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so.*

%files -n libbpp-phyl-omics-devel
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%dir %{_prefix}/lib/cmake/
%dir %{_prefix}/lib/cmake/bpp-phyl-omics
%{_prefix}/%{_lib}/lib*.so
%{_prefix}/%{_lib}/lib*.a
%{_prefix}/lib/cmake/bpp-phyl-omics/bpp-phyl-omics*.cmake
%{_prefix}/include/*

%changelog
* Wed May 10 2017 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.3.0-1
- Added distance matrix output filter.
- Added clock and reparametrization options in model fitting.
* Fri Sep 26 2014 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.2.0-1
- Added model fitting and parameter estimations.
* Thu Mar 07 2013 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.1.0-1
- Initial release. Contains tools to build phylogenies along a genome alignment (distance methods only for now).
* Tue Nov 06 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.3-1
- First draft of the spec file.

