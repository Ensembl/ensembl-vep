# Ensembl VEP — Local Installation Notes

## Environment

- **VEP version**: 115.2
- **Platform**: macOS (Apple Silicon / arm64)
- **Perl**: 5.42 (via Homebrew)
- **Install date**: 2026-05-27

---

## Installation Steps

### 1. Clone the repository

```zsh
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
git checkout release/115
```

### 2. Install Perl dependencies

```zsh
# Install cpanminus if not present
brew install cpanminus

# Install DBI
cpanm DBI

# Install DBD::mysql (requires mysql-client)
brew install mysql-client

# Build DBD::mysql manually with correct library paths
cpanm --look DBD::mysql
# Inside the source dir:
LIBRARY_PATH="/opt/homebrew/opt/openssl@3/lib:/opt/homebrew/opt/zstd/lib:/opt/homebrew/opt/mysql-client/lib" \
perl Makefile.PL --mysql_config=/opt/homebrew/opt/mysql-client/bin/mysql_config
make && make install
```

### 3. Download caches

```zsh
# GRCh37 (hg19)
perl INSTALL.pl --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh37 --DESTDIR ~/.vep --NO_UPDATE

# GRCh38 (hg38)
perl INSTALL.pl --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh38 --DESTDIR ~/.vep --NO_UPDATE
```

> **Note**: Caches are ~15GB each. For reliable downloads, run in a persistent terminal using `nohup`:
> ```zsh
> nohup perl INSTALL.pl --AUTO cf --SPECIES homo_sapiens --ASSEMBLY GRCh37 --DESTDIR ~/.vep --NO_UPDATE > ~/vep_install.log 2>&1 &
> tail -f ~/vep_install.log
> ```

### 4. Add VEP to your shell environment

Add to `~/.zshrc`:

```zsh
# Ensembl VEP
export PATH="$HOME/ensembl-vep:$PATH"
export PERL5LIB="$HOME/.vep:$HOME/ensembl-vep/modules:$PERL5LIB"
```

Then reload: `source ~/.zshrc`

---

## Cache Locations

| Assembly | Path |
|----------|------|
| GRCh37   | `~/.vep/homo_sapiens/115_GRCh37` |
| GRCh38   | `~/.vep/homo_sapiens/115_GRCh38` |

---

## Usage

### Basic annotation (GRCh37)

```zsh
vep \
  --input_file input.vcf \
  --output_file output.txt \
  --cache \
  --assembly GRCh37 \
  --species human \
  --dir_cache ~/.vep \
  --force_overwrite
```

### Basic annotation (GRCh38)

```zsh
vep \
  --input_file input.vcf \
  --output_file output.txt \
  --cache \
  --assembly GRCh38 \
  --species human \
  --dir_cache ~/.vep \
  --force_overwrite
```

### Common flags

| Flag | Description |
|------|-------------|
| `--vcf` | Output in VCF format |
| `--everything` | Enable all annotations (SIFT, PolyPhen, AF, etc.) |
| `--sift b` | SIFT predictions (score + label) |
| `--polyphen b` | PolyPhen predictions |
| `--af_gnomade` | gnomAD exome allele frequencies |
| `--pick` | One consequence per variant (most severe) |
| `--fork 4` | Parallel processing with 4 CPUs |
| `--stats_file stats.html` | Generate summary statistics HTML |

### Filter results

```zsh
filter_vep \
  --input_file output.txt \
  --filter "IMPACT is MODERATE or IMPACT is HIGH"
```
