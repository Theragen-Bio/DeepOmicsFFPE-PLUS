# Changelog

## 0.3.1-rc2 (2026-04-01)

### 변경 사항

#### `--ref-version` 인자 제거 및 `--ref-fasta` 로 대체
- `doffpe_client.py`: `--ref-version` (`-r`) 인자 제거, `--ref-fasta` 인자로 레퍼런스 FASTA 경로를 직접 받도록 변경
- `ffpe_utils/01_ext_bam.py`: 동일하게 `--ref-version` 제거, `--ref-fasta` 사용으로 변경

#### 버전 업데이트 (`0.3.0-rc8` → `0.3.1-rc2`)
- `doffpe/__init__.py`: `__version__` 업데이트
- `pyproject.toml`: `version` 업데이트
- `setup.cfg`: `version` 업데이트

#### 잔존 stale 텍스트 수정
- `doffpe_client.py`: `--help` epilog 예제에서 `-r hg38` 제거 → `--ref-fasta /path/to/hg38.fa` 로 변경
- `ffpe_utils/01_ext_bam.py`: `--help` epilog 예제에서 `-r hg38` 제거 → `--ref-fasta /path/to/hg38.fa` 로 변경
- `ffpe_utils/01_ext_bam.py`: `--ref-fasta` 도움말에서 `"must match --ref-version"` 문구 제거

#### 버그 수정
- `ffpe_utils/01_ext_bam.py`: `finally` 블록에서 `after_bam_time`, `end_time` 미초기화로 인한 `NameError` 방지
  - `try` 블록 진입 전에 `after_bam_time = end_time = start_time` 초기화 추가

> 서버 쪽 변경사항은 `dofp_pipeline/CHANGES.md` 참조
